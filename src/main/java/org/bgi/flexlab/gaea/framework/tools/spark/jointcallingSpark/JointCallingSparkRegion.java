package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.Locus;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.hdfs.util.ByteArray;
import org.apache.spark.Accumulator;
import org.apache.spark.Partition;
import org.apache.spark.SparkConf;
import org.apache.spark.TaskContext;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.*;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.serializer.JavaSerializer;
import org.apache.spark.serializer.Serializer;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.annotator.interval.Genome;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcallingprepare.WindowsBasedTestPartitioner;
import org.bgi.flexlab.gaea.util.Utils;
import org.jcodings.util.Hash;
import org.seqdoop.hadoop_bam.VariantContextWithHeader;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Array;
import scala.Tuple2;

import java.io.*;
import java.net.URI;
import java.rmi.Remote;
import java.util.*;
import java.util.zip.GZIPInputStream;

import static org.apache.spark.api.java.StorageLevels.MEMORY_ONLY;

public class JointCallingSparkRegion {
    public final static String INPUT_ORDER = "input.name.order";
    public final static String INPUT_LIST = "input.gvcf.list";// added by gc
    public final static String Window_File = "window.file.path";// window file path
    public static Logger logger = LoggerFactory.getLogger(JointCallingSparkRegion.class);
    public static LinkedHashMap<String, Integer> chrIndex=new LinkedHashMap<>();
    //    public static ArrayList<BufferedReader> sampleReaders=new ArrayList<>();
//    public static ArrayList<BufferedReader> sampleReadersForMR2=new ArrayList<>();//HashMap可以深拷贝，不担心互相受影响
    public static GenomeLocation parseRegionFromString(String targetRegion) {
//        GenomeLocation gloc;
        String ele="";
        ArrayList<String> eles=new ArrayList<>();
        for(int i=0;i<targetRegion.length();i++){
            if(Character.isLetterOrDigit(targetRegion.charAt(i))){
                ele+=targetRegion.charAt(i);
            }else{
                eles.add(ele);
                ele="";
                if(eles.size()==3){
                    break;
                }
            }
        }
        if(eles.size()==2){
            if(ele!=""){
                eles.add(ele);
            }
        }
        if(eles.size()!=3){
            logger.error("fail to parse targetRegion, please check the options");
        }
        String chr=eles.get(0);
        if(chrIndex.containsKey(chr)){
            return new GenomeLocation(chr,chrIndex.get(chr),Integer.parseInt(eles.get(1)),Integer.parseInt(eles.get(2)));
        }else if(chrIndex.containsValue(Integer.parseInt(chr))){
            for(Map.Entry<String,Integer> kv:chrIndex.entrySet()){
                if(kv.getValue()==Integer.parseInt(chr)){
                    return new GenomeLocation(kv.getKey(),kv.getValue(),Integer.parseInt(eles.get(1)),Integer.parseInt(eles.get(2)));
                }
            }
            logger.error("code error");
        }else{
            logger.error("no such chromosome name/id:"+chr+", please check the options");
        }
        logger.error("fail to parse targetRegion, please check the options");
        return null;
    }
    public static void main(String[] args) throws IOException, ClassNotFoundException {
        String HEADER_DEFAULT_PATH = "vcfheader";
        String MERGER_HEADER_INFO = "vcfheaderinfo";
        ArrayList<ArrayList<String>> multiMapSampleNames=new ArrayList<>();
        LinkedHashMap<String,VCFFileReader> sampleReader=new LinkedHashMap<>();

        JointCallingSparkOptions options = new JointCallingSparkOptions();

        String vVcfPath;
        LinkedHashMap<String,Integer> sampleIndex=new LinkedHashMap<>();
        Map<String, String> pathSample = new HashMap<>();


        //OutputStream win_out = null;
        VCFHeader header = null;
        LinkedHashMap<Integer, String> contigs = null;
        SparkConf conf = new SparkConf().setAppName("JointCallingSparkRegion");

        conf.set("spark.serializer", "org.apache.spark.serializer.JavaSerializer");
        conf.set("spark.rdd.compress","true");
//        conf.set("spark.kryo.registrator","org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark.JointCallingKryoRegistrator");
        JavaSparkContext sc = new JavaSparkContext(conf);   //打开spark环境
        Configuration hadoopConf=null;
        hadoopConf=sc.hadoopConfiguration();
        options.parse(args);//获得header等信息
        String outputDir=options.getOutDir();
        //生成虚拟vcf header
        File tmpDir=new File(options.getTmpOut());
        if(!tmpDir.exists()){
            tmpDir.mkdirs();
        }
        Path vf=new Path(options.getOutDir()+"/virtual.vcf");
        FileSystem vfs=vf.getFileSystem(hadoopConf);
        BufferedWriter vfOut = new BufferedWriter(new OutputStreamWriter(vfs.create(vf)));
        String fileList=options.getInputList();
        if(fileList.startsWith("file://")){
            fileList=fileList.substring(7);
        }
        BufferedReader gvcfList=new BufferedReader(new FileReader(fileList));

        String gvcfFilePath;
        if((gvcfFilePath=gvcfList.readLine())!=null){
            BufferedReader gvcfReader=null;
            if(gvcfFilePath.startsWith("file://")) {
                gvcfFilePath=gvcfFilePath.substring(7);
                if (gvcfFilePath.endsWith(".gz")) {
                    gvcfReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfFilePath))));
                } else {
                    gvcfReader = new BufferedReader(new FileReader(gvcfFilePath));
                }
            }else{
                FileSystem fs = FileSystem.get(URI.create(gvcfFilePath),hadoopConf);
                FSDataInputStream fsr = fs.open(new Path(gvcfFilePath));
                if(gvcfFilePath.endsWith(".gz")){
                    gvcfReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(fsr)));
                }else {
                    gvcfReader = new BufferedReader(new InputStreamReader(fsr));
                }
            }
            String tmpline;
            while ((tmpline = gvcfReader.readLine()) != null) {
                if (tmpline.startsWith("#CHROM")) {
                    String writeLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVirtualSample\n";
                    vfOut.write(writeLine);
                    break;
                } else {
                    vfOut.write(tmpline);
                    vfOut.write("\n");
                }
            }
            gvcfReader.close();
            vfOut.close();
        }
        //使用虚拟vcf文件，提取contig信息
        SeekableStream in=WrapSeekable.openPath(vfs, vf);
        header = VCFHeaderReader.readHeaderFrom(in);
        in.close();
        if (header == null)
            throw new RuntimeException("header is null !!!");
        contigs = new LinkedHashMap<>();
        for (VCFContigHeaderLine line : header.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
            chrIndex.put(line.getID(),line.getContigIndex());
        }

        //创建窗口文件
        String win_out_file = options.getOutDir() + "/windows.bed";
        // Path raw_win_file=new Path(options.getWinFile());//get windows file path from
        // command
        Path winPath=new Path(win_out_file);
        FileSystem fs=winPath.getFileSystem(hadoopConf);
        if (!fs.exists(winPath)) {
            fs.createNewFile(winPath);
        }
        BufferedWriter win_out = new BufferedWriter(new OutputStreamWriter(fs.create(winPath))); // output stream ready for write
        int window_size = options.getWindowsSize();
        if(options.getTargetRegion()!=null){
            GenomeLocation targetRegion=parseRegionFromString(options.getTargetRegion());
            int start=targetRegion.getStart();
            int end=-1;
            while (true) {
                end = start + window_size - 1;
                if (end > targetRegion.getEnd()) {
                    end = targetRegion.getEnd();
                }
                String write_line = targetRegion.getContig() + "\t" + start + "\t" + end + "\n";
                win_out.write(write_line);
                start += window_size;
                if (start > targetRegion.getEnd()) {
                    break;
                }
            }
        }else {
            for (Map.Entry<Integer, String> entry : contigs.entrySet()) {
                // System.out.println(entry.getKey()+"\t"+entry.getValue());
                String chr = entry.getValue();
                int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
                int start = 1;
                int end = -1;
                while (true) {
                    end = start + window_size - 1;
                    if (end > contigLength) {
                        end = contigLength;
                    }
                    String write_line = chr + "\t" + start + "\t" + end + "\n";
                    win_out.write(write_line);
                    start += window_size;
                    if (start > contigLength) {
                        break;
                    }
                }
            }
        }
        win_out.close();
        JavaRDD<String> windows=sc.textFile(win_out_file,options.getMapperNumber());

        Configuration hadoop_conf=new Configuration(hadoopConf);
        MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
        logger.warn("before get Header");
        ArrayList<Path> pathList=new ArrayList<>();
        for(String a:options.getInputStringList()){
            pathList.add(new Path(a));
        }
        if (options.getVcfHeaderFile() != null) {
        } else {
            multiVcfHeader.headersConfig(pathList, options.getVCFHeaderOutput(), hadoop_conf);
        }
        JavaRDD<String> gvcfSamples = sc.textFile(options.getInputList(), options.getMapperNumber());
        Path outDir=new Path(options.getOutDir()+"/index");
        FileSystem outDirFs=outDir.getFileSystem(hadoopConf);
        if(!outDirFs.exists(outDir)){
            outDirFs.mkdirs(outDir);
            JavaPairRDD<Integer,String> offsets=gvcfSamples.flatMapToPair(new buildIndex(options));
            offsets.groupByKey().foreach(new printIndex(options));
        }
        windows.foreachPartition(new GenotypeGVCFs(options));
        sc.stop();
        return;
//        HashMap<String,String> confMap=new HashMap<>();
//        confMap.put(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
//        confMap.put(MERGER_HEADER_INFO, options.getOutput()+"/"+MERGER_HEADER_INFO);
//        ArrayList<String> inputListArr = new ArrayList<>();
//        JavaRDD<String> indexInfo=sc.textFile(hadoop_conf.get(MERGER_HEADER_INFO),1).persist(MEMORY_ONLY);
//
//
//        for(String info:indexInfo.collect()){
//            String[] eles = info.split("\t");
//            if (eles.length != 3) {
//                logger.error("vcfheaderinfo file format error");
//            }
//            String name;
//            if (eles[2].endsWith(",")) {
//                name = eles[2].substring(0, eles[2].length() - 1);
//            } else {
//                name = eles[2];
//            }
//            sampleIndex.put(name,Integer.valueOf(eles[1]));
//            pathSample.put(eles[0], name);
//        }
//        indexInfo.unpersist();
//        gvcfList.close();
//
//        boolean virtualDone=false;
//        int bufReaderSize=1024;
//        //当样本量大时，这个bufferReader就容易造成OOM错误
//        long startMem = Runtime.getRuntime().freeMemory();
//        int iter2=0;
//        long endMem = Runtime.getRuntime().freeMemory();
//        double useMemory = (startMem - endMem)/1024.00/1024.00;
//        System.out.println("使用掉的内存："+ useMemory + "M");
//        System.out.println("内存总量："+Runtime.getRuntime().totalMemory()/1024/1024 + "M");
////        System.exit(1);
//        gvcfList=new BufferedReader(new FileReader(fileList));
//        String tmpLine;
//        while ((tmpLine = gvcfList.readLine()) != null) {
//            String[] eles = tmpLine.split("/");
//            String rName = eles[eles.length - 1];
//            String sampleName = pathSample.get(rName);
//            inputListArr.add(sampleName);
//        }
//        gvcfList.close();
//        String[] inputListStringArr = new String[inputListArr.size()];
//        for (int i = 0; i < inputListArr.size(); i++) {
//            inputListStringArr[i] = inputListArr.get(i);
//        }
//        hadoop_conf.set(INPUT_ORDER, Utils.join(",", inputListStringArr));
//        confMap.put(INPUT_ORDER, Utils.join(",", inputListStringArr));
//
//        inputListArr.clear();
//        hadoop_conf.set(INPUT_LIST, options.getInputList());
//        confMap.put(INPUT_LIST, options.getInputList());
//        // create windows bed file based on window size
//        hadoop_conf.set(Window_File, options.getTmpOut() + "/windows.bed");
//        confMap.put(Window_File, options.getTmpOut() + "/windows.bed");
//        vVcfPath=options.getTmpOut()+"/virtual.vcf";
//        logger.warn("after get Header");
//
//        Integer chrInt=0;
//        Integer startPos=1;
//        Integer endPos=50000000;
//        //在这之前必须要获取每条染色体长度，来判断终止位置是否到了染色体末尾
//        JavaRDD<String> gvcfSamples = sc.textFile(options.getInputList(), options.getMapperNumber());
//        int[] index=new int[options.getMapperNumber()];
//        for(int i=0;i<options.getMapperNumber();i++){
//            index[i]=i;
//        }
//        for(List<String> samples:gvcfSamples.collectPartitions(index)){
//            ArrayList<String> indexSamples=new ArrayList<>();
//            for(String sample:samples){
//                String[] eles = sample.split("/");
//                String rName = eles[eles.length - 1];
//                String sampleName = pathSample.get(rName);
//                indexSamples.add(sampleName);
//            }
//            multiMapSampleNames.add(indexSamples);
//        }
//        int iter=0;
//        Broadcast<DriverBC> dBC=sc.broadcast(new DriverBC(outputDir,multiMapSampleNames,chrIndex,options, vVcfPath, sampleIndex, pathSample));
//        while(true){
//            if(chrInt>=chrIndex.size()){
//                break;
//            }
//            String chr=contigs.get(chrInt);
//            GenomeLocation region=new GenomeLocation(chr,chrInt,startPos,endPos);
//            if(options.getTargetRegion()!=null){
//                GenomeLocation targetRegion=parseRegionFromString(options.getTargetRegion());
//                region=targetRegion;
//            }
//            JavaPairRDD<GenomeLocation,Integer> variantsRegion = gvcfSamples.flatMapToPair(new ProcessVariantLocus(region,dBC));
//            JavaPairRDD<GenomeLocation,Integer> partitionedRegion=variantsRegion.partitionBy(new GenomeLocPartitioner(options.getMapperNumber(),region));
//            JavaRDD<GenomeLocation> mergedRegion=partitionedRegion.keys().mapPartitions(new MergeRegion(region)).repartition(1);
//            String outputBP=options.getOutDir()+"/bps."+String.valueOf(iter);
//            mergedRegion.saveAsTextFile(outputBP);
//            Path bpDir=new Path(outputBP);
//            FileSystem fs=bpDir.getFileSystem(hadoopConf);
//            int codeInspectorIter=0;
//            if(fs.isDirectory(bpDir)){
//                RemoteIterator<LocatedFileStatus> files=fs.listFiles(bpDir,true);
//                while(files.hasNext()){
//                    LocatedFileStatus file=files.next();
//                    String[] names=file.getPath().toString().split("/");
//                    String lastName=names[names.length-1];
//                    if(lastName.startsWith("part")){
//                        codeInspectorIter++;
//                    }
//                }
//            }
//            if(codeInspectorIter>1){
//                logger.error("code error, only one BPs file should be generated");
//                System.exit(1);
//            }
//            String mergedAllBPs="";
//            if(fs.isDirectory(bpDir)){
//                RemoteIterator<LocatedFileStatus> files=fs.listFiles(bpDir,true);
//                while(files.hasNext()){
//                    LocatedFileStatus file=files.next();
//                    String[] names=file.getPath().toString().split("/");
//                    String lastName=names[names.length-1];
//                    if(lastName.startsWith("part")){
//                        mergedAllBPs=file.getPath().toString();
//                    }
//                }
//            }
//            JavaRDD<VariantContext> combinedVariants2= gvcfSamples.mapPartitions(new CombineVariants(region,mergedAllBPs,confMap,dBC));
//            JavaPairRDD<Integer, VariantContext> combinedVariants=combinedVariants2.mapToPair(new PairFunction<VariantContext, Integer, VariantContext>() {
//                @Override public Tuple2<Integer, VariantContext> call(VariantContext variantContext) throws Exception {
//                    return new Tuple2<>(variantContext.getStart(),variantContext);
//                }
//            });
//            JavaPairRDD<Integer, VariantContext> partitionedCombineVariants=combinedVariants.repartitionAndSortWithinPartitions(new GenomeLocPartitioner2(options.getReducerNumber(),region));
//            JavaRDD<String> genotypedVariants=partitionedCombineVariants.values().mapPartitions(new GenotypeVariants(region,args,mergedAllBPs,confMap,dBC));
//            genotypedVariants.saveAsTextFile(options.getOutput()+"/genotypeVCF."+String.valueOf(iter));
//            iter++;
//            if(options.getTargetRegion()!=null){
//                chrInt=chrIndex.size()+1;
//            }
//        }
//        if(tmpDir.isDirectory()){
//            String[] files=tmpDir.list();
//            for(String file:files){
//                if(file.endsWith(".lastVC")){
//                    (new File(file)).delete();
//                }
//            }
//        }
//        sc.stop();
    }
}
