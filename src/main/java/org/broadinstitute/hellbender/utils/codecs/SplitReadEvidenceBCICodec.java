package org.broadinstitute.hellbender.utils.codecs;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class SplitReadEvidenceBCICodec extends AbstractBCICodec<SplitReadEvidence> {
    private boolean versionChecked = false;
    private static final String SR_BCI_FILE_EXTENSION = ".sr.bci";

    @Override
    public SplitReadEvidence decode( final Reader<SplitReadEvidence> reader ) throws IOException {
        if ( !versionChecked ) {
            if ( !SplitReadEvidence.BCI_VERSION.equals(reader.getVersion()) ) {
                throw new UserException("baf.bci file has wrong version: expected " +
                        SplitReadEvidence.BCI_VERSION + " but found " + reader.getVersion());
            }
            versionChecked = true;
        }
        final DataInputStream dis = reader.getStream();
        final String sample = reader.getSampleNames().get(dis.readInt());
        final String contig = reader.getDictionary().getSequence(dis.readInt()).getSequenceName();
        final int position = dis.readInt();
        final int count = dis.readInt();
        final boolean strand = dis.readBoolean();
        return new SplitReadEvidence(sample, contig, position, count, strand);
    }

    @Override
    public Class<SplitReadEvidence> getFeatureType() { return SplitReadEvidence.class; }

    @Override
    public boolean canDecode( final String path ) { return path.endsWith(SR_BCI_FILE_EXTENSION); }

    public void encode( final SplitReadEvidence srEvidence,
                        final Writer<SplitReadEvidence> writer ) throws IOException {
        final DataOutputStream dos = writer.getStream();
        dos.writeInt(writer.getSampleIndex(srEvidence.getSample()));
        dos.writeInt(writer.getContigIndex(srEvidence.getContig()));
        dos.writeInt(srEvidence.getStart());
        dos.writeInt(srEvidence.getCount());
        dos.writeBoolean(srEvidence.getStrand());
    }

    public String getVersion() { return SplitReadEvidence.BCI_VERSION; }
}
