package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import com.esotericsoftware.kryo.Kryo;
import de.javakaffee.kryoserializers.protobuf.ProtobufSerializer;
import org.apache.spark.serializer.KryoRegistrator;

public class JointCallingKryoRegistrator implements KryoRegistrator {
    @Override public void registerClasses(Kryo kryo) {
        kryo.register(htsjdk.variant.variantcontext.FastGenotype.class,new ProtobufSerializer(),10100);
//        kryo.register(htsjdk.variant.variantcontext.Allele.class,new ProtobufSerializer(),10101);
//        kryo.register(htsjdk.variant.variantcontext.GenotypesContext.class,new ProtobufSerializer(),10102);
    }
}
