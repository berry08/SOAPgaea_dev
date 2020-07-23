package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.VoidFunction;
import scala.Tuple2;

import java.io.BufferedWriter;
import java.io.OutputStreamWriter;
import java.util.Iterator;

public class printIndex implements VoidFunction<Tuple2<Integer, Iterable<String>>> {
    public JointCallingSparkOptions options;
    public printIndex(JointCallingSparkOptions options) {
        this.options=options;
    }

    @Override public void call(Tuple2<Integer, Iterable<String>> integerIterableTuple2) throws Exception {
        Path outPath=new Path(options.getOutDir()+"/index/"+integerIterableTuple2._1);
        Configuration conf=new Configuration();
        FileSystem fs=outPath.getFileSystem(conf);
        BufferedWriter out=new BufferedWriter(new OutputStreamWriter(fs.create(outPath)));
        Iterator<String> values=integerIterableTuple2._2.iterator();
        while(values.hasNext()){
            String curString=values.next();
//            String[] eles=curString.split("\t");
//            String sampleName=eles[0];
//            Long offset=Long.parseLong(eles[1]);
            out.write(curString+"\n");
        }
        out.close();
    }
}
