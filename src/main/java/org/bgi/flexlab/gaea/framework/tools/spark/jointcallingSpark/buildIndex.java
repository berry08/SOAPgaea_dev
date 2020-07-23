package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.api.java.function.VoidFunction;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import scala.Tuple2;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class buildIndex implements PairFlatMapFunction<String, Integer, String>{
        public JointCallingSparkOptions options;
        public buildIndex(JointCallingSparkOptions options) {
            this.options=options;
        }

        @Override public Iterator<Tuple2<Integer, String>> call(String s) throws Exception {
            Configuration conf=new Configuration();
            Path virualHeaderPath=new Path(options.getOutDir()+"/virtual.vcf");
            SeekableStream in2 = WrapSeekable.openPath(virualHeaderPath.getFileSystem(conf), virualHeaderPath);
            VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in2);
            Set<VCFHeaderLine> gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
            VCFHeaderVersion version=null;
            HashMap<String,Integer> chrIndex=new HashMap<>();
            HashMap<Integer,String> contigs=new HashMap<>();
            for (VCFContigHeaderLine line : virtualHeader.getContigLines()) {
                contigs.put(line.getContigIndex(), line.getID());
                chrIndex.put(line.getID(),line.getContigIndex());
            }
            long[] accumulateLength=new long[chrIndex.size()];
            long totalLength=0;
            for(int ii=0;ii<chrIndex.size();ii++){
                long chrLength=virtualHeader.getSequenceDictionary().getSequence(ii).getSequenceLength();
                totalLength+=chrLength;
                accumulateLength[ii]=totalLength-chrLength;
            }
            for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
                if (VCFHeaderVersion.isFormatString(line.getKey())) {
                    version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                    break;
                }
            }
            ArrayList<Tuple2<Integer,String>> returnValue=new ArrayList<>();
            String sample=s;
            Path samplePath=new Path(sample);
            String sampleName=null;
            if(sample.startsWith("file://")) {
                sample=sample.substring(7);
            }
            BufferedReader reader=null;
            if(sample.endsWith(".gz")) {
                reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(sample))));
            }else {
                reader = new BufferedReader(new InputStreamReader(new FileInputStream(sample)));
            }
            String line;
            VCFCodec codec=null;
            long curOffset=0;
            while(true) {
                line = reader.readLine();
                if (line.startsWith("##")) {
                    continue;
                } else if (line.startsWith("#")) {
                    String[] eles=line.split("\t");
                    sampleName=eles[eles.length-1];
                    codec=new VCFCodec();
                    ArrayList<String> curSamples=new ArrayList<>();
                    curSamples.add(eles[eles.length-1]);
                    virtualHeader=new VCFHeader(gvcfHeaderMetaInfo,curSamples);
                    codec.setVCFHeader(virtualHeader,version);
                    break;
                }
            }



            long skipLength=1024*500;
            if(reader.skip(skipLength)<skipLength){
                return returnValue.iterator();
            }
            curOffset+=skipLength;
            int lastPosM=0;
            while((line=reader.readLine())!=null){
                curOffset+=line.length()+1;
                if((line=reader.readLine())!=null){
                    VariantContext v=codec.decode(line);
                    long curPos=accumulateLength[chrIndex.get(v.getContig())]+v.getStart();
                    int curPosM= (int) (curPos/1000000);
                    if(curPosM!=lastPosM){
                        String value=sampleName+"\t"+curOffset;
                        for(int i=lastPosM+1;i<=curPosM;i++) {
                            returnValue.add(new Tuple2<>(i, value));
                        }
                    }
                    lastPosM=curPosM;
                    curOffset+=line.length()+1;
                }
                if(reader.skip(skipLength)<skipLength){
                    break;
                }
                curOffset+=skipLength;
            }
            System.out.println(sampleName+" index build done");
            return returnValue.iterator();
        }

}
