package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.*;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Header;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

/**
 * Prints SV evidence records. Can be used with -L to retrieve records on a set of intervals. Supports streaming input
 * from GCS buckets.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Coordinate-sorted and indexed evidence file URI
 *     </li>
 *     <li>
 *         Reference sequence dictionary
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Coordinate-sorted evidence file, automatically indexed if ending with ".gz"
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk PrintSVEvidence \
 *       --evidence-file gs://my-bucket/batch_name.SR.txt.gz \
 *       -L intervals.bed \
 *       --sequence-dictionary ref.dict \
 *       -O local.SR.txt.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Prints SV evidence records to a file",
        oneLineSummary = "Prints SV evidence records",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
@DocumentedFeature
public final class PrintSVEvidence extends FeatureWalker<Feature> {

    public static final String EVIDENCE_FILE_NAME = "evidence-file";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";

    @Argument(
            doc = "Input file URI with extension '"
                    + SplitReadEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + DiscordantPairEvidenceCodec.FORMAT_SUFFIX + "', '"
                    + BafEvidenceCodec.FORMAT_SUFFIX + "', or '"
                    + DepthEvidenceCodec.FORMAT_SUFFIX + "' (may be gzipped). "
                    + "Can also handle bci rather than txt files.",
            fullName = EVIDENCE_FILE_NAME
    )
    private GATKPath inputFilePath;

    @Argument(
            doc = "Output file with an evidence extension matching the input. Will be indexed if it has a " +
                    "block-compressed extension (e.g. '.gz' or '.bci').",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    @Argument(doc = "List of sample names", fullName = "sample-names", optional = true)
    private List<String> sampleNames = new ArrayList<>();

    @SuppressWarnings("rawtypes")
    private FeatureOutputStream foStream;
    @SuppressWarnings("rawtypes")
    private Writer bciWriter;
    private FeatureCodec<? extends Feature, ?> featureCodec;
    private Class<? extends Feature> evidenceClass;

    private static final List<FeatureCodec<? extends Feature, ?>> outputCodecs = new ArrayList<>(8);
    static {
        outputCodecs.add(new BafEvidenceCodec());
        outputCodecs.add(new DepthEvidenceCodec());
        outputCodecs.add(new DiscordantPairEvidenceCodec());
        outputCodecs.add(new SplitReadEvidenceCodec());
        outputCodecs.add(new BafEvidenceBCICodec());
        outputCodecs.add(new DepthEvidenceBCICodec());
        outputCodecs.add(new DiscordantPairEvidenceBCICodec());
        outputCodecs.add(new SplitReadEvidenceBCICodec());
    }

    @Override
    protected boolean isAcceptableFeatureType(final Class<? extends Feature> featureType) {
        return featureType.equals(BafEvidence.class) || featureType.equals(DepthEvidence.class)
                || featureType.equals(DiscordantPairEvidence.class) || featureType.equals(SplitReadEvidence.class);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        return inputFilePath;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        featureCodec = FeatureManager.getCodecForFile(inputFilePath.toPath());
        evidenceClass = featureCodec.getFeatureType();
        initializeOutput();
    }

    private static FeatureCodec<? extends Feature, ?> findOutputCodec( final GATKPath outputFilePath ) {
        final String outputFileName = outputFilePath.toString();
        for ( final FeatureCodec<? extends Feature, ?> codec : outputCodecs ) {
            if ( codec.canDecode(outputFileName) ) {
                return codec;
            }
        }
        throw new UserException("no codec found for path " + outputFileName);
    }

    private void initializeOutput() {
        final FeatureCodec<? extends Feature, ?> outputCodec = findOutputCodec(outputFilePath);
        final Class<? extends Feature> outputClass = outputCodec.getFeatureType();
        if ( evidenceClass != outputClass ) {
            throw new UserException("input file contains " + evidenceClass.getSimpleName() +
                    " but output file would be expected to contain " + outputClass.getSimpleName());
        }

        if ( outputCodec instanceof AbstractBCICodec ) {
            final AbstractBCICodec<? extends Feature> bciCodec =
                    (AbstractBCICodec<? extends Feature>)outputCodec;
            final Object headerObj = getDrivingFeaturesHeader();
            final Header header;
            if ( headerObj instanceof Header ) {
                header = (Header)headerObj;
            } else {
                header = new Header(outputClass.getSimpleName(), bciCodec.getVersion(),
                                   getBestAvailableSequenceDictionary(), sampleNames);
            }
            bciWriter = new Writer<>(outputFilePath, header, bciCodec::encode, compressionLevel);
        } else {
            if ( evidenceClass.equals(DiscordantPairEvidence.class) ) {
                initializeFOStream(new DiscordantPairEvidenceCodec(), DiscordantPairEvidenceCodec::encode);
            } else if ( evidenceClass.equals(SplitReadEvidence.class) ) {
                initializeFOStream(new SplitReadEvidenceCodec(), SplitReadEvidenceCodec::encode);
            } else if ( evidenceClass.equals(BafEvidence.class) ) {
                initializeFOStream(new BafEvidenceCodec(), BafEvidenceCodec::encode);
            } else if ( evidenceClass.equals(DepthEvidence.class) ) {
                initializeFOStream(new DepthEvidenceCodec(), DepthEvidenceCodec::encode);
                writeDepthEvidenceHeader();
            } else {
                throw new UserException.BadInput("Unsupported evidence type: " + evidenceClass.getSimpleName());
            }
        }
    }

    private <F extends Feature> void initializeFOStream( final FeatureCodec<F, ?> codec,
                                                         final Function<F, String> encoder ) {
        foStream = new FeatureOutputStream<>(outputFilePath, codec, encoder,
                getBestAvailableSequenceDictionary(), compressionLevel);
    }

    private void writeDepthEvidenceHeader() {
        final Object header = getDrivingFeaturesHeader();
        if ( header instanceof Header ) {
            final StringBuilder sb = new StringBuilder("#Chr\tStart\tEnd");
            final Header hdr = (Header)header;
            for ( final String sampleName : hdr.getSampleNames() ) {
                sb.append('\t').append(sampleName);
            }
            foStream.writeHeader(sb.toString());
        } else if (header instanceof String) {
            foStream.writeHeader((String) header);
        } else {
            throw new IllegalArgumentException("Expected header object of type String");
        }
    }

    @Override
    @SuppressWarnings("unchecked")
    public void apply(final Feature feature,
                      final ReadsContext readsContext,
                      final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        if ( foStream != null ) {
            foStream.add(feature);
        }
        if ( bciWriter != null ) {
            bciWriter.write(feature);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if ( foStream != null ) {
            foStream.close();
        }
        if ( bciWriter != null ) {
            bciWriter.close();
        }
        return null;
    }
}
