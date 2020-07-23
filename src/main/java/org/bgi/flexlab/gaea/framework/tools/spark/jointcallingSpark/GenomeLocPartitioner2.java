package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.apache.spark.Partitioner;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

public class GenomeLocPartitioner2 extends Partitioner {
    public int partitions;
    public GenomeLocation gloc;
    public int partitionRange;
    public GenomeLocPartitioner2(int parts,GenomeLocation region){
        this.partitions=parts;
        this.gloc=region;
        partitionRange=(gloc.getEnd()-gloc.getStart())/partitions;
    }
    @Override
    public int numPartitions() {
        return partitions;
    }
    /*
    todo
     */
    @Override
    public int getPartition(Object key) {
        Integer glocKey=(Integer) key;
//        if(glocKey.getStart()<1100000){
//            return 0;
//        }else {
//            return 1;
//        }
        if(glocKey-gloc.getStart()<0){
            return 0;
        }else if(glocKey>gloc.getEnd()){
            return partitions-1;
        }else {
            return (glocKey - gloc.getStart()) / partitionRange;
        }
    }
}
