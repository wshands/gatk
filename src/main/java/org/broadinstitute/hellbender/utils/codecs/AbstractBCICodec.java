package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.samtools.util.LocationAware;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Writer;

import java.io.IOException;
import java.io.InputStream;

import static htsjdk.tribble.FeatureCodecHeader.EMPTY_HEADER;

public abstract class AbstractBCICodec <T extends Feature> implements FeatureCodec<T, Reader<T>> {

    @Override
    public Feature decodeLoc( final Reader<T> reader ) throws IOException {
        return decode(reader);
    }

    @Override
    public FeatureCodecHeader readHeader( final Reader<T> reader ) throws IOException {
        return EMPTY_HEADER;
    }

    @Override
    public Reader<T> makeSourceFromStream( final InputStream is ) {
        throw new GATKException("wasn't expecting to execute this code path");
    }

    @Override
    public LocationAware makeIndexableSourceFromStream( final InputStream is ) {
        throw new GATKException("wasn't expecting to execute this code path");
    }

    @Override
    public boolean isDone( final Reader<T> reader ) { return !reader.hasNext(); }

    @Override
    public void close( final Reader<T> reader ) { reader.close(); }

    abstract public void encode( T feature, Writer<T> writer ) throws IOException;

    abstract public String getVersion();
}
