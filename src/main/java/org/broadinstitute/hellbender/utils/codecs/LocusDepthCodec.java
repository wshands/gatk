package org.broadinstitute.hellbender.utils.codecs;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.walkers.sv.PairedEndAndSplitReadEvidenceCollection.LocusDepth;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.IOException;

public class LocusDepthCodec extends AbstractBCICodec<LocusDepth> {
    private boolean versionChecked = false;
    private static final String LD_BCI_FILE_EXTENSION = ".ld.bci";

    @Override
    public LocusDepth decode( final Reader<LocusDepth> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !LocusDepth.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("bci file has wrong version: expected " +
                        LocusDepth.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        return new LocusDepth(reader.getDictionary(), reader.getStream());
    }

    @Override
    public Class<LocusDepth> getFeatureType() { return LocusDepth.class; }

    @Override
    public boolean canDecode( final String path ) { return path.endsWith(LD_BCI_FILE_EXTENSION); }

    public void encode( final LocusDepth locusDepth, final Writer<LocusDepth> writer )
            throws IOException {
        locusDepth.write(writer.getStream());
    }

    public String getVersion() { return LocusDepth.BCI_VERSION; }
}
