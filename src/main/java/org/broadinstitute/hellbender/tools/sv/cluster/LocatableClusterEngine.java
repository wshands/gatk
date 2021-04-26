package org.broadinstitute.hellbender.tools.sv.cluster;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
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

    @VisibleForTesting
    abstract boolean clusterTogether(final T a, final T b);
    abstract protected int getMaxClusterableStartingPosition(final T item);

    public final List<T> getOutput() {
        flushClusters();
        final List<T> output = new ArrayList<>(outputBuffer);
        outputBuffer.clear();
        return output;
    }

    public final boolean isEmpty() {
        return idToClusterMap.isEmpty() && outputBuffer.isEmpty();
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

    private final int registerItem(final T item) {
        final int itemId = nextItemId++;
        idToItemMap.put(itemId, item);
        return itemId;
    }

    private final int getMaxClusterableStartingPositionByIds(final Collection<Integer> itemIds) {
        Utils.nonNull(itemIds);
        Utils.nonEmpty(itemIds);
        return getMaxClusterableStartingPosition(itemIds.stream().map(this::getItem).collect(Collectors.toList()));
    }

    @VisibleForTesting
    protected final int getMaxClusterableStartingPosition(final Collection<T> items) {
        final List<String> contigA = items.stream().map(T::getContigA).distinct().collect(Collectors.toList());
        if (contigA.size() > 1) {
            throw new IllegalArgumentException("Items start on multiple contigs");
        }
        return items.stream().mapToInt(item -> getMaxClusterableStartingPosition(item)).max().getAsInt();
    }

    @VisibleForTesting
    protected final Function<Collection<T>, T> getCollapser() {
        return collapser;
    }

    /**
     * Add a new {@param <T>} to the current clusters and determine which are complete
     * @param itemId
     * @return the IDs for clusters that are complete and ready for processing
     */
    private final List<Integer> cluster(final Integer itemId) {
        final T item = getItem(itemId);
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
            final List<Integer> clusterItems = cluster.getItemIds();
            if (item.getPositionA() > cluster.getMaxClusterableStart()) {
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
        idToClusterMap.put(nextClusterId++, new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems));
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
        idToClusterMap.put(nextClusterId++, new Cluster(getMaxClusterableStartingPosition(getItem(item)), newClusters));
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
        final Cluster newCluster = new Cluster(getMaxClusterableStartingPositionByIds(newClusterItems), newClusterItems);

        //Do not add duplicates
        if (!idToClusterMap.entrySet().contains(newCluster)) {
            idToClusterMap.put(nextClusterId++, newCluster);
        }
    }

    private final Cluster getCluster(final int id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Cluster ID " + id + " does not exist.");
        return idToClusterMap.get(id);
    }

    private final T getItem(final int id) {
        Utils.validateArg(idToItemMap.containsKey(id), "Item ID " + id + " does not exist.");
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
        final List<Integer> clusterItems = cluster.getItemIds();
        clusterItems.add(itemId);
        final T item = getItem(itemId);
        final int itemClusterableStartPosition = getMaxClusterableStartingPosition(item);
        cluster.setMaxClusterableStart(Math.max(cluster.getMaxClusterableStart(), itemClusterableStartPosition));
    }

    private static final class Cluster {
        private int maxClusterableStart;
        private final List<Integer> itemIds;

        public Cluster(final int maxClusterableStart, final List<Integer> itemIds) {
            Utils.nonNull(itemIds);
            this.maxClusterableStart = maxClusterableStart;
            this.itemIds = itemIds;
        }

        public int getMaxClusterableStart() {
            return maxClusterableStart;
        }

        public void setMaxClusterableStart(final int position) {
            maxClusterableStart = position;
        }

        public List<Integer> getItemIds() {
            return itemIds;
        }
    }

}
