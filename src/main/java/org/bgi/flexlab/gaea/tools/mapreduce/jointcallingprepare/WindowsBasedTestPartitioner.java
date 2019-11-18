package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare;

import org.apache.hadoop.mapreduce.Partitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;

public class WindowsBasedTestPartitioner<T> extends Partitioner<WindowsBasedWritable, T> {

	@Override
	public int getPartition(WindowsBasedWritable key, T v, int numPartitioner) {
		//chr1 1M-2M
		Integer pos=key.getPosition().get();
		int band=(int)(1000000*1.1/numPartitioner);
		return (int)((pos-1000000) / band);
	}
}