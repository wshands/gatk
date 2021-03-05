package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public abstract class LocatableClusterEngine<T extends SVLocatable> {

    public enum CLUSTERING_TYPE {
        SINGLE_LINKAGE,
        MAX_CLIQUE
    }

    protected final SAMSequenceDictionary dictionary;
    private final Function<Collection<T>, T> collapser; // Flattens clusters into a single representative item for output
    private Map<Integer, Cluster> idToClusterMap; // Active clusters
    private final List<T> outputBuffer;
    protected final CLUSTERING_TYPE clusteringType;
    private String currentContig;
    private final Map<Integer, T> idToItemMap;
    private int nextItemId;
    private int nextClusterId;


    public LocatableClusterEngine(final SAMSequenceDictionary dictionary,
                                  final CLUSTERING_TYPE clusteringType,
                                  final Function<Collection<T>, T> collapser) {
        this.dictionary = Utils.nonNull(dictionary);
        this.clusteringType = clusteringType;
        this.collapser = Utils.nonNull(collapser);
        idToClusterMap = new HashMap<>();
        outputBuffer = new ArrayList<>();
        currentContig = null;
        idToItemMap = new HashMap<>();
        nextItemId = 0;
        nextClusterId = 0;
    }

    abstract protected boolean clusterTogether(final T a, final T b);
    abstract protected SimpleInterval getClusteringInterval(final T item);

    protected SimpleInterval getClusteringInterval(final T item, final SimpleInterval clusterMinStartInterval) {
        final SimpleInterval itemInterval = getClusteringInterval(item);
        if (clusterMinStartInterval == null) {
            return itemInterval;
        } else {
            return restrictIntervalsForClustering(itemInterval, clusterMinStartInterval);
        }
    }

    protected final SimpleInterval restrictIntervalsForClustering(final SimpleInterval a, final SimpleInterval b) {
        Utils.validateArg(a.getContig().equals(b.getContig()), "Intervals are on different contigs");
        final int start = Math.min(a.getStart(), b.getStart());
        final int end = Math.max(a.getEnd(), b.getEnd());
        return new SimpleInterval(a.getContig(), start, end);
    }

    protected final SimpleInterval getClusteringInterval(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        final List<T> items = itemIds.stream().map(this::getItem).collect(Collectors.toList());
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        final List<SimpleInterval> clusteringIntervals = items.stream().map(item -> getClusteringInterval(item, null)).collect(Collectors.toList());
        final int minStart = clusteringIntervals.stream().mapToInt(SimpleInterval::getStart).min().getAsInt();
        final int maxEnd = clusteringIntervals.stream().mapToInt(SimpleInterval::getEnd).max().getAsInt();
        return new SimpleInterval(contigA.get(0), minStart, maxEnd);
    }

    public final List<T> getOutput() {
        flushClusters();
        final List<T> output = new ArrayList<>(outputBuffer);
        outputBuffer.clear();
        return output;
    }

    public final boolean isEmpty() {
        return outputBuffer.isEmpty();
    }

    private final int registerItem(final T item) {
        final int itemId = nextItemId++;
        idToItemMap.put(itemId, item);
        return itemId;
    }

    public final void add(final T item) {
        // Start a new cluster if on a new contig
        if (!item.getContigA().equals(currentContig)) {
            flushClusters();
            currentContig = item.getContigA();
            seedCluster(registerItem(item));
            return;
        }
        final int itemId = registerItem(item);
        final List<Integer> clusterIdsToProcess = cluster(itemId);
        processFinalizedClusters(clusterIdsToProcess);
    }

    public final String getCurrentContig() {
        return currentContig;
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId
     * @return the IDs for clusters that are complete and ready for processing
     */
    private final List<Integer> cluster(final Integer itemId) {
        final T item = getItem(itemId);
        if (item.getType().equals(StructuralVariantType.DEL) && item.getPositionA() == 11212835) {
            int x = 0;
        }
        // Get list of item IDs from active clusters that cluster with this item
        final Set<Integer> linkedItems = idToClusterMap.values().stream().map(Cluster::getItemIds)
                .flatMap(List::stream)
                .distinct()
                .filter(other -> !other.equals(itemId) && clusterTogether(item, getItem(other)))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        // Find clusters to which this item belongs, and which active clusters we're definitely done with
        final List<Integer> clusterIdsToProcess = new ArrayList<>();
        final List<Integer> clustersToAugment = new ArrayList<>();
        final Set<List<Integer>> clustersToSeedWith = new HashSet<>();    // Use set to prevent creating duplicate clusters
        for (final Map.Entry<Integer, Cluster> entry : idToClusterMap.entrySet()) {
            final Integer clusterIndex = entry.getKey();
            final Cluster cluster = entry.getValue();
            final SimpleInterval clusterInterval = cluster.getInterval();
            final List<Integer> clusterItems = cluster.getItemIds();
            if (getClusteringInterval(item, null).getStart() > clusterInterval.getEnd()) {
                clusterIdsToProcess.add(clusterIndex);  //this cluster is complete -- process it when we're done
            } else {
                if (clusteringType.equals(CLUSTERING_TYPE.MAX_CLIQUE)) {
                    final List<Integer> linkedClusterItems = clusterItems.stream().filter(linkedItems::contains).collect(Collectors.toList());
                    final int numLinkedItems = linkedClusterItems.size();
                    if (numLinkedItems == clusterItems.size()) {
                        clustersToAugment.add(clusterIndex);
                    } else if (numLinkedItems > 0) {
                        clustersToSeedWith.add(linkedClusterItems);
                    }
                } else if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
                    final boolean matchesCluster = clusterItems.stream().anyMatch(linkedItems::contains);
                    if (matchesCluster) {
                        clustersToAugment.add(clusterIndex);
                    }
                } else {
                    throw new IllegalArgumentException("Clustering algorithm for type " + clusteringType.name() + " not implemented");
                }
            }
        }

        // Create new clusters from subsets (max-clique only)
        if (!clustersToSeedWith.isEmpty()) {
            final List<Set<Integer>> augmentedClusterItemLists = clustersToAugment.stream()
                    .map(idToClusterMap::get)
                    .map(Cluster::getItemIds)
                    .map(HashSet::new)
                    .collect(Collectors.toList());
            final List<Set<Integer>> triggeredClusterItemSets = new ArrayList<>(clustersToSeedWith.size() + augmentedClusterItemLists.size());
            triggeredClusterItemSets.addAll(clustersToSeedWith.stream().map(HashSet::new).collect(Collectors.toList()));
            triggeredClusterItemSets.addAll(augmentedClusterItemLists);
            triggeredClusterItemSets.sort(Comparator.comparingInt(Set::size));
            for (int i = 0; i < clustersToSeedWith.size(); i++) {
                final Set<Integer> seedItems = triggeredClusterItemSets.get(i);
                // Check that this cluster is not a sub-cluster of any of the others being created
                boolean isSubset = false;
                for (int j = i + 1; j < triggeredClusterItemSets.size(); j++) {
                    if (isSubsetOf(seedItems, triggeredClusterItemSets.get(j))) {
                        isSubset = true;
                        break;
                    }
                }
                if (!isSubset) {
                    seedWithExistingCluster(itemId, seedItems);
                }
            }
        }

        // Add to or merge existing clusters
        if (clusteringType.equals(CLUSTERING_TYPE.SINGLE_LINKAGE)) {
            if (!clustersToAugment.isEmpty()) {
                combineClusters(clustersToAugment, itemId);
            }
        } else {
            for (final Integer clusterId : clustersToAugment) {
                addToCluster(clusterId, itemId);
            }
        }

        // If there weren't any matches, create a new singleton cluster
        if (clustersToAugment.isEmpty() && clustersToSeedWith.isEmpty()) {
            seedCluster(itemId);
        }
        return clusterIdsToProcess;
    }

    private final void combineClusters(final Collection<Integer> clusterIds, final Integer itemId) {
        final List<Cluster> clusters = clusterIds.stream().map(this::getCluster).collect(Collectors.toList());
        clusterIds.stream().forEach(idToClusterMap::remove);
        final List<Integer> clusterItems = clusters.stream()
                .map(Cluster::getItemIds)
                .flatMap(List::stream)
                .distinct()
                .collect(Collectors.toList());
        final List<Integer> newClusterItems = new ArrayList<>(clusterItems.size() + 1);
        newClusterItems.addAll(clusterItems);
        newClusterItems.add(itemId);
        idToClusterMap.put(nextClusterId++, new Cluster(getClusteringInterval(newClusterItems), newClusterItems));
    }

    private final void processCluster(final int clusterIndex) {
        final Cluster cluster = getCluster(clusterIndex);
        idToClusterMap.remove(clusterIndex);
        outputBuffer.add(collapser.apply(cluster.getItemIds().stream().map(idToItemMap::get).collect(Collectors.toList())));
    }

    private final void processFinalizedClusters(final List<Integer> clusterIdsToProcess) {
        for (final Integer clusterId : clusterIdsToProcess) {
            processCluster(clusterId);
        }
    }

    private final boolean isSubsetOf(final Set<Integer> smallerSet, final Set<Integer> largerSet) {
        if (largerSet.size() < smallerSet.size()) {
            return false;
        }
        int n = 0;
        for (final Integer item : smallerSet) {
            if (largerSet.contains(item)) {
                if (++n == smallerSet.size()) {
                    return true;
                }
            }
        }
        return false;
    }

    private final void flushClusters() {
        final List<Integer> clustersToFlush = new ArrayList<>(idToClusterMap.keySet());
        for (final Integer clusterId : clustersToFlush) {
            processCluster(clusterId);
        }
        idToItemMap.clear();
        nextItemId = 0;
        nextClusterId = 0;
    }

    private final void seedCluster(final Integer item) {
        final List<Integer> newClusters = new ArrayList<>(1);
        newClusters.add(item);
        idToClusterMap.put(nextClusterId++, new Cluster(getClusteringInterval(getItem(item), null), newClusters));
    }

    /**
     * Create a new cluster
     * @param item
     * @param seedItems
     */
    private final void seedWithExistingCluster(final Integer item, final Collection<Integer> seedItems) {
        final List<Integer> newClusterItems = new ArrayList<>(1 + seedItems.size());
        newClusterItems.addAll(seedItems);
        newClusterItems.add(item);
        final Cluster newCluster = new Cluster(getClusteringInterval(getItem(item), null), newClusterItems);

        //Do not add duplicates
        if (!idToClusterMap.entrySet().contains(newCluster)) {
            idToClusterMap.put(nextClusterId++, newCluster);
        }
    }

    protected final Cluster getCluster(final int id) {
        if (!idToClusterMap.containsKey(id)) {
            throw new IllegalArgumentException("Specified cluster ID " + id + " does not exist.");
        }
        return idToClusterMap.get(id);
    }

    protected final T getItem(final int id) {
        if (!idToItemMap.containsKey(id)) {
            throw new IllegalArgumentException("Specified item ID " + id + " does not exist.");
        }
        return idToItemMap.get(id);
    }

    /**
     * Add the item specified by {@param itemId} to the cluster specified by {@param clusterIndex}
     * and expand the clustering interval
     * @param clusterId
     * @param itemId
     */
    private final void addToCluster(final int clusterId, final Integer itemId) {
        final Cluster cluster = getCluster(clusterId);
        final T item = getItem(itemId);
        final SimpleInterval clusterInterval = cluster.getInterval();
        final List<Integer> clusterItems = cluster.getItemIds();
        clusterItems.add(itemId);
        final SimpleInterval clusteringStartInterval = getClusteringInterval(item, clusterInterval);
        if (clusteringStartInterval.getStart() != clusterInterval.getStart() || clusteringStartInterval.getEnd() != clusterInterval.getEnd()) {
            idToClusterMap.put(clusterId, new Cluster(clusteringStartInterval, clusterItems));
        }
    }

    private static final class Cluster {
        private final SimpleInterval interval;
        private final List<Integer> itemIds;

        public Cluster(final SimpleInterval interval, final List<Integer> itemIds) {
            Utils.nonNull(interval);
            Utils.nonNull(itemIds);
            this.interval = interval;
            this.itemIds = itemIds;
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public List<Integer> getItemIds() {
            return itemIds;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o.getClass().equals(Cluster.class))) return false;
            Cluster cluster = (Cluster) o;
            return Objects.equals(interval, cluster.interval) && Objects.equals(itemIds, cluster.itemIds);
        }

        @Override
        public int hashCode() {
            return Objects.hash(interval, itemIds);
        }
    }

}
