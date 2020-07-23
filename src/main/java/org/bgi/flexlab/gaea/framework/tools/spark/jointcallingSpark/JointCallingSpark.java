package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.Locus;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
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
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalWriter;
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

public class JointCallingSpark {
    public final static String INPUT_ORDER = "input.name.order";
    public final static String INPUT_LIST = "input.gvcf.list";// added by gc
    public final static String Window_File = "window.file.path";// window file path
    public static Logger logger = LoggerFactory.getLogger(JointCallingSpark.class);
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
        SparkConf conf = new SparkConf().setAppName("JointCallingSpark");

        conf.set("spark.serializer", "org.apache.spark.serializer.JavaSerializer");
        conf.set("spark.rdd.compress","true");
//        conf.set("spark.kryo.registrator","org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark.JointCallingKryoRegistrator");
        JavaSparkContext sc = new JavaSparkContext(conf);   //打开spark环境
        Configuration hadoopConf=null;
        hadoopConf=sc.hadoopConfiguration();
        options.parse(args);//获得header等信息
        String outputDir=options.getOutDir();
        //生成虚拟vcf header
        File tmpDir=new File(options.getOutDir());
        if(!tmpDir.exists()){
            tmpDir.mkdirs();
        }


        Configuration hadoop_conf=new Configuration(hadoopConf);
        MultipleVCFHeaderForJointCalling multiVcfHeader = new MultipleVCFHeaderForJointCalling();
        logger.warn("before get Header");
//        ArrayList<Path> pathList=new ArrayList<>();
        VCFHeader mergeHeader=null;
//        for(String a:options.getInputStringList()){
//            pathList.add(new Path(a));
//        }
        File headerDir=new File(options.getOutDir()+"/headers");
        if(!headerDir.exists()){
            headerDir.mkdirs();
        }else{
            if(!headerDir.isDirectory()){
                headerDir.delete();
                headerDir.mkdirs();
            }
        }

        File vcfheaderFile=new File(options.getOutDir()+"/vcfheader");
        File vcfheaderInfoFile=new File(options.getOutDir()+"/vcfheaderinfo");
        if(!vcfheaderFile.exists() || !vcfheaderInfoFile.exists() || !headerDir.exists()) {
            JavaRDD<String> gvcfSamples = sc.textFile(options.getInputList(), options.getMapperNumber());
            gvcfSamples.mapPartitionsWithIndex(new ProcessHeader(options), true).collect();
        }
        //create vcfheader and vcfheaderinfo by small files in directory "headers"
        BufferedWriter vcfInfoWriter=new BufferedWriter(new FileWriter(options.getOutDir()+"/"+MERGER_HEADER_INFO));
        int inputIndex=0;
        LinkedHashSet<VCFHeader> headers=new LinkedHashSet<>();
        TreeSet<String> sampleNames=new TreeSet<>();
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
        hadoop_conf.set(MERGER_HEADER_INFO, options.getOutDir()+"/"+MERGER_HEADER_INFO);
        HashMap<String,String> sampleNamePath=new HashMap<>();
        HashMap<String,String> sampleNameRPath=new HashMap<>();
        for(int i=0;i<options.getMapperNumber();i++){
            File smallHeader=new File(options.getOutDir()+"/headers/vcfheader"+i);
            File smallHeaderInfo=new File(options.getOutDir()+"/headers/vcfPathName"+i);
            if(!smallHeader.exists()){
                logger.error("file not exists,"+smallHeader.getAbsolutePath());
                System.exit(1);
            }
            if(!smallHeaderInfo.exists()){
                logger.error("file not exists,"+smallHeaderInfo.getAbsolutePath());
                System.exit(1);
            }

            BufferedReader smallHeaderInfoReader=new BufferedReader(new FileReader(smallHeaderInfo.getAbsolutePath()));
            String headerInfoLine;
            while((headerInfoLine=smallHeaderInfoReader.readLine())!=null){
                String[] eles=headerInfoLine.split("\t");
                sampleNamePath.put(eles[2],eles[0]);
                sampleNameRPath.put(eles[2],eles[1]);
//                vcfHeaderWriter.write(eles[0]+"\t"+inputIndex+"\t"+eles[2]+"\n");
//                sampleIndex.put(eles[2],inputIndex);
//                pathSample.put(eles[0], eles[2]);
//                inputIndex++;
            }
            smallHeaderInfoReader.close();
            VCFLocalLoader vcfLL=new VCFLocalLoader(smallHeader.getAbsolutePath());
            headers.add(vcfLL.getHeader());
            sampleNames.addAll(vcfLL.getHeader().getSampleNamesInOrder());
            if(inputIndex>=1000){
                VCFHeader vcfHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
                headers.clear();
                headers.add(vcfHeader);
            }
        }
        if (headers.size() >= 1) {
            mergeHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
            headers.clear();
        }
        inputIndex=0;
        File rawInput=new File(options.getInputList());
        String sortedInputGvcfList=options.getOutDir()+"/sorted."+rawInput.getName();
        BufferedWriter sortedInputWriter=new BufferedWriter(new FileWriter(sortedInputGvcfList));
        for(String sampleName:mergeHeader.getSampleNamesInOrder()){
            if(!sampleNamePath.containsKey(sampleName)){
                logger.error("code error");
                System.exit(1);
            }
            sortedInputWriter.write(sampleNamePath.get(sampleName)+"\n");
            vcfInfoWriter.write(sampleNameRPath.get(sampleName)+"\t"+inputIndex+"\t"+sampleName+"\n");
            sampleIndex.put(sampleName,inputIndex);
            pathSample.put(sampleNameRPath.get(sampleName), sampleName);
            inputIndex++;
        }
        sortedInputWriter.close();
        vcfInfoWriter.close();
        VCFLocalWriter mergedVCFHeaderWriter=new VCFLocalWriter(options.getOutDir()+"/"+HEADER_DEFAULT_PATH,false,false);
        mergedVCFHeaderWriter.writeHeader(mergeHeader);
        mergedVCFHeaderWriter.close();
        sortedInputGvcfList="file://"+sortedInputGvcfList;
//        options.setInputList(sortedInputGvcfList);
        JavaRDD<String> sortedGvcfSamples = sc.textFile(sortedInputGvcfList, options.getMapperNumber());
        FileWriter virtualVCF=new FileWriter(new File(options.getOutDir()+"/virtual.vcf"));
        String fileList=sortedInputGvcfList;
        if(fileList.startsWith("file://")){
            fileList=fileList.substring(7);
        }
        BufferedReader headerReader=new BufferedReader(new FileReader(options.getOutDir()+"/vcfheader"));

        String tmpline;
        while ((tmpline = headerReader.readLine()) != null) {
            if (tmpline.startsWith("#CHROM")) {
                String writeLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVirtualSample\n";
                virtualVCF.write(writeLine);
                break;
            } else {
                virtualVCF.write(tmpline);
                virtualVCF.write("\n");
            }
        }
        virtualVCF.close();
        //使用虚拟vcf文件，提取contig信息
        SeekableStream in=new SeekableFileStream(new File(options.getOutDir()+"/virtual.vcf"));
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
        if (win_out_file.startsWith("file://")) {
            win_out_file = win_out_file.substring(7);
        }
        File winOutFile = new File(win_out_file);
        if (!winOutFile.exists()) {
            winOutFile.createNewFile();
        }
        FileWriter win_out = new FileWriter(winOutFile); // output stream ready for write
        int window_size = options.getWindowsSize();
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
        win_out.close();
//        System.exit(1);
        HashMap<String,String> confMap=new HashMap<>();
//        confMap.put(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
        confMap.put(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getOutDir() + "/"+HEADER_DEFAULT_PATH);
//        confMap.put(MERGER_HEADER_INFO, options.getOutput()+"/"+MERGER_HEADER_INFO);
        confMap.put(MERGER_HEADER_INFO, options.getOutDir()+"/"+MERGER_HEADER_INFO);
        ArrayList<String> inputListArr = new ArrayList<>();

//        Path vcfHeaderInfo = new Path(hadoop_conf.get(MERGER_HEADER_INFO));
//        FileSystem vhiFs = vcfHeaderInfo.getFileSystem(hadoop_conf);
//        BufferedReader indexReader = new BufferedReader()(new InputStreamReader(vhiFs.open(vcfHeaderInfo)));

        int bufReaderSize=1024;
        //当样本量大时，这个bufferReader就容易造成OOM错误
        long startMem = Runtime.getRuntime().freeMemory();
        BufferedReader gvcfList=new BufferedReader(new FileReader(fileList));
        String tmpLine;
        while ((tmpLine = gvcfList.readLine()) != null) {
            String[] eles = tmpLine.split("/");
            String rName = eles[eles.length - 1];
            String sampleName = pathSample.get(rName);
            inputListArr.add(sampleName);
        }
        gvcfList.close();
        String[] inputListStringArr = new String[inputListArr.size()];
        for (int i = 0; i < inputListArr.size(); i++) {
            inputListStringArr[i] = inputListArr.get(i);
        }
        hadoop_conf.set(INPUT_ORDER, Utils.join(",", inputListStringArr));
        confMap.put(INPUT_ORDER, Utils.join(",", inputListStringArr));

        inputListArr.clear();
        hadoop_conf.set(INPUT_LIST, sortedInputGvcfList);
        confMap.put(INPUT_LIST, sortedInputGvcfList);
        // create windows bed file based on window size
        hadoop_conf.set(Window_File, options.getOutDir() + "/windows.bed");
        confMap.put(Window_File, options.getOutDir() + "/windows.bed");
        vVcfPath=options.getOutDir()+"/virtual.vcf";
        logger.warn("after get Header");

        //获取变异断点，不分区域处理，因为最后所有位置要集中分区，不然最后结果无法做到排序整合。这里内存做一个限制，因为当样本量巨大（百万）时，假设每个样本变异数量为3M,那么总内存大约为 (4+4+4)*3M*1M=36T,这一步如果我们把内存限制在10G，那么程序可以做个智能判断，针对样品数量进行个评估，如果少于1000个，不管是否时WGS，都不拆分区域；如果超过1000，则检查区域大小（参数指定，默认为WGS），如果在允许范围内则不分区域，否则按照窗口区域进行处理
        //2020/03/30 目前的思路时按照超大样本WGS的逻辑来写，即分区域来处理
        //如果将样本分区，对每个分区按区域提取位置，并合并，当把所有区域提取出来后，再跟其它分区一起重分区排序，这样还是会存在巨大样本量时潜在的内存问题，还是要回到分区域思路上，一次处理50M（这个大小可以保证单个分区变异位置合并后重分区的总大小在100G以内，大部分集群都可以满足），生成这50M的变异区间，然后在下游做变异合并

        //获得每个样本的reader，这是老代码，适用于普通文件，不适用于hdfs系统
//        gvcfList=new BufferedReader()(new FileReader(fileList));
//        String gvcfPath;
//        while((gvcfPath=gvcfList.readLine())!=null){
//            VCFFileReader tmp=new VCFFileReader(new File(gvcfPath));
//            sampleReader.put(gvcfPath,tmp);
//        }

        //由于VCFFileReader不支持Path类型，所以处理不了hdfs文件，在hadoop_bam包未找到对应的类，因此在各个partition中使用query方法无法实现。其实从hdfs的特点也可以得知，很难实现一个像query的函数，这个时候索引就显得意义不大了，所以改成用BufferedReader的方式，由于在executor中循环顺序处理各个区域，而BufferReader又不能seek，将会导致效率很低，所以只能放在driver中，当样本量非常大时再去控制buffer大小（默认是8K，1M的样本也才8G，暂时能接受，处理好尾巴问题就行）
        //在上面生成virtual.vcf时，同时生成所有文件的BufferReader

        //分区域处理，每个区域执行一次SPARK
        //当用户未指定处理区域时，这里分区大小改为50M，不再从磁盘文件读取，直接在内存里生成
        //在这之前必须要获取每条染色体长度，来判断终止位置是否到了染色体末尾

        int[] index=new int[options.getMapperNumber()];
        for(int i=0;i<options.getMapperNumber();i++){
            index[i]=i;
        }
        for(List<String> samples:sortedGvcfSamples.collectPartitions(index)){
            ArrayList<String> indexSamples=new ArrayList<>();
            for(String sample:samples){
                String[] eles = sample.split("/");
                String rName = eles[eles.length - 1];
                String sampleName = pathSample.get(rName);
                indexSamples.add(sampleName);
            }
            multiMapSampleNames.add(indexSamples);
        }
        int iter=0;

//        Broadcast<ArrayList<BufferedReader>> sampleReadersForMR2BC=sc.broadcast(sampleReadersForMR2);
        HashMap<Integer,Long> accumulateLength=new HashMap<>();
        long totalLength=0;
        for(int ii=0;ii<chrIndex.size();ii++){
            long chrLength=header.getSequenceDictionary().getSequence(ii).getSequenceLength();
            totalLength+=chrLength;
            accumulateLength.put(ii,totalLength-chrLength);
        }
        SeekableStream in2=new SeekableFileStream(new File(options.getOutDir()+"/virtual.vcf"));
        VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in2);
        Set<VCFHeaderLine> gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        VCFHeaderVersion version=null;
        for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
            if (VCFHeaderVersion.isFormatString(line.getKey())) {
                version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                break;
            }
        }

        Broadcast<DriverBC> dBC=sc.broadcast(new DriverBC(options.getOutDir(),multiMapSampleNames,chrIndex,options, vVcfPath, sampleIndex, pathSample,accumulateLength,virtualHeader,version));
        Long step=50000000L-1;
        Long longStart=1L;
        Long longEnd=longStart+step;
        GenomeLongRegion region=new GenomeLongRegion(longStart,longEnd);
        Long cycleEnd=totalLength;
        if(options.getTargetRegion()!=null){
            GenomeLocation targetRegion=parseRegionFromString(options.getTargetRegion());
            region.setStart(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getStart());
            if(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd()-accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getStart()>step){
                region.setEnd(region.getStart()+step);
            }else{
                region.setEnd(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd());
            }
            if(accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd()<cycleEnd){
                cycleEnd=accumulateLength.get(chrIndex.get(targetRegion.getContig()))+targetRegion.getEnd();
            }
        }

        while(true){
            //结束条件1：基因组区域处理结束
            if(region.getStart()>cycleEnd){
                break;
            }
            ArrayList<GenomeLocation> regions=new ArrayList<>();
            long startLoc=0;
            long endLoc=0;
            ArrayList<Integer> spanChrs=new ArrayList<>();
            int startChr=-1;
            int endChr=-1;
            for(int ii=0;ii<accumulateLength.size();ii++){
                if(accumulateLength.get(ii)>region.getStart()){
                    if(startChr==-1) {
                        startChr=ii-1;
                        startLoc = region.getStart() - accumulateLength.get(ii - 1);
                    }
                }
                if(accumulateLength.get(ii)>=region.getEnd()){
                    if(endChr==-1) {
                        endChr=ii-1;
                        endLoc = region.getEnd() - accumulateLength.get(ii - 1);
                        break;
                    }
                }
            }
            if(endChr==startChr){
                regions.add(new GenomeLocation(contigs.get(startChr),(int)startLoc,(int)endLoc));
            }else{
                for(int ii=startChr;ii<=endChr;ii++){
                    if(ii==startChr){
                        regions.add(new GenomeLocation(contigs.get(ii),(int)startLoc,header.getSequenceDictionary().getSequence(ii).getSequenceLength()));
                    }else if(ii==endChr){
                        regions.add(new GenomeLocation(contigs.get(ii),1,(int)endLoc));
                    }else{
                        regions.add(new GenomeLocation(contigs.get(ii),1,header.getSequenceDictionary().getSequence(ii).getSequenceLength()));
                    }
                }
            }

//            String chr=contigs.get(chrInt);

            //partition必须依赖PairRDD，所以下面必须转换成PairRDD
            //提取变异位置
            JavaPairRDD<GenomeLongRegion,Integer> variantsRegion = sortedGvcfSamples.flatMapToPair(new ProcessVariantLocus(region,regions,dBC));
            //分区
            String outputBP=options.getOutDir()+"/bps."+String.valueOf(iter);
            File bpDir=new File(outputBP);
            if(!bpDir.exists() || !bpDir.isDirectory()) {

                JavaPairRDD<GenomeLongRegion, Integer> partitionedRegion = variantsRegion.partitionBy(new GenomeLocPartitioner(options.getMapperNumber(), region));
//            System.out.println(variantsRegion.toDebugString());
                //每个分区内部整合和排序，最后合并成一个，这里暂时不写到磁盘
//            System.out.println(partitionedRegion.getNumPartitions()+"\t"+partitionedRegion.partitioner());
//            System.out.println(posList.getNumPartitions()+"\t"+posList.partitioner());
                JavaRDD<GenomeLongRegion> mergedRegion = partitionedRegion.keys().mapPartitions(new MergeRegion(region)).sortBy(x -> x.getStart(), true, 1);
                //也可以同时分区和排序
                //JavaPairRDD<GenomeLocation,Integer> sortedRegion=partitionedRegion.repartitionAndSortWithinPartitions(new GenomeLocPartitioner(5,region),new GenomeLocComp());
                //分区内去重

                //生成潜在变异区间文件
                mergedRegion.saveAsTextFile("file://" + outputBP);
            }
//            System.exit(1);
//            File bpDir=new File(outputBP);
            int codeInspectorIter=0;
            if(bpDir.isDirectory()){
                String[] files=bpDir.list();
                for(String file:files){
                    if(file.startsWith("part")){
                        codeInspectorIter++;
                    }
                }
            }
            if(codeInspectorIter>1){
                logger.error("code error, only one BPs file should be generated");
                System.exit(1);
            }
            String mergedAllBPs="";
            if(bpDir.isDirectory()){
                String[] files=bpDir.list();
                for(String file:files){
                    if(file.startsWith("part")){
                        mergedAllBPs=bpDir.getAbsolutePath()+"/"+file;
                    }
                }
            }
            Path pBPpath=new Path(mergedAllBPs);
            BufferedReader bp_reader=new BufferedReader(new FileReader(mergedAllBPs));
            int bpIter=0;
            String bpRegion=bp_reader.readLine();
            bpIter++;
            ArrayList<Long> bpDis=new ArrayList<>();
            if(bpRegion==null){
                bpIter++;
                region.setStart(region.getEnd()+1);
                region.setEnd(region.getStart()+step);
                continue;
            }else{
                while((bpRegion=bp_reader.readLine())!=null){
                    bpIter++;
                }
            }
            bp_reader.close();
            int variantsNumInEachReducer=bpIter/options.getReducerNumber();
            bp_reader=new BufferedReader(new FileReader(mergedAllBPs));
            ArrayList<Long> bpPartition=new ArrayList<>();
            bpIter=0;
            String lastBpRegion=null;
            while((bpRegion=bp_reader.readLine())!=null){
                bpIter++;
                if(bpIter%variantsNumInEachReducer==0 && bpPartition.size()<options.getReducerNumber()-1){
                    String[] eles=bpRegion.split("\t");
                    bpPartition.add(Long.parseLong(eles[0]));
                }
                lastBpRegion=bpRegion;
            }
            bp_reader.close();
            String[] bpEles=lastBpRegion.split("\t");
            bpPartition.add(Long.parseLong(bpEles[0])+1);
            if(bpPartition.size()!=options.getReducerNumber()){
                System.out.println("code error: bp partition error");
                System.exit(1);
            }
            //做第二步MR
            //考虑到内存问题（在mapreduce中是shuffle），还是将样本分批，并且先做combineGVCFs，这样会减少内存中流动数据的大小，然后再做genotypeGVCFs
//            JavaPairRDD<Integer, VariantContext> combinedVariants= sortedGvcfSamples.mapPartitionsToPair(new CombineVariants(region,mergedAllBPs,confMap,dBC));
//            JavaRDD<VariantContext> combinedVariants2= sortedGvcfSamples.mapPartitions(new CombineVariants(region,mergedAllBPs,confMap,dBC));
            //hdfs路径
//            String combineOut=options.getOutDir()+"/combine."+String.valueOf(iter);
            //本地路径
            String combineOutLocal=options.getOutDir()+"/combine."+String.valueOf(iter);
            //生成潜在变异区间文件
//            Path combineOutPath=new Path(combineOut);
//            FileSystem combineFs=combineOutPath.getFileSystem(hadoopConf);
//            if(!combineFs.exists(combineOutPath)){
//                combineFs.mkdirs(combineOutPath);
//            }else{
//                if(!combineFs.isDirectory(combineOutPath)){
//                    combineFs.delete(combineOutPath);
//                }
//            }
            File combineOutLocalFile=new File(combineOutLocal);
            if(!combineOutLocalFile.exists()){
                combineOutLocalFile.mkdirs();
                sortedGvcfSamples.foreachPartition(new CombineVariants(region,regions,mergedAllBPs,confMap,dBC,iter,bpPartition));
            }else{
                if(!combineOutLocalFile.isDirectory()){
                    combineOutLocalFile.delete();
                    combineOutLocalFile.mkdirs();
                    sortedGvcfSamples.foreachPartition(new CombineVariants(region,regions,mergedAllBPs,confMap,dBC,iter,bpPartition));
                }
            }

//            JavaPairRDD<Integer, Integer> combinedVariants= sortedGvcfSamples.mapPartitionsToPair(new ExtractVariants(region,mergedAllBPs,confMap,dBC,iter));
            //这一步涉及到大量数据的网络传输，是限速步骤之一
//            String genotypeOut=options.getOutDir()+"/genotype."+String.valueOf(iter);
            String genotypeOut=options.getOutDir()+"/genotype."+String.valueOf(iter);
            //生成潜在变异区间文件
            File genotypeOutPath=new File(genotypeOut);
            if(!genotypeOutPath.exists()){
                genotypeOutPath.mkdirs();
            }else{
                if(!genotypeOutPath.isDirectory()){
                    genotypeOutPath.delete();
                }
            }
            JavaRDD<String> noneGvcfSamples=sortedGvcfSamples.repartition(options.getReducerNumber());
            JavaRDD<String> variantsNum=noneGvcfSamples.mapPartitionsWithIndex(new ParseCombineAndGenotypeGVCFs(region,regions,args,mergedAllBPs,confMap,dBC,iter,bpPartition),true);
            variantsNum.collect();
//            JavaPairRDD<Integer, VariantContext> partitionedCombineVariants=combinedVariants.repartitionAndSortWithinPartitions(new GenomeLocPartitioner2(options.getReducerNumber(),region));
//            JavaRDD<String> genotypedVariants=partitionedCombineVariants.values().mapPartitions(new GenotypeVariants(region,args,mergedAllBPs,confMap,dBC));
//            JavaRDD<String> genotypedVariants=partitionedCombineVariants.values().mapPartitionsWithIndex(new GenotypeVariantsNoMerge(region,args,mergedAllBPs,confMap,dBC),true);

//            partitionedCombineVariants.saveAsTextFile(options.getOutput()+"/combineVCF."+String.valueOf(iter));
            //下一行的机制有点奇怪，当combinedVariants保存到磁盘时，会导致genotypedVariants输出结果为空，待以后查证依据
//            combinedVariants.saveAsTextFile(options.getOutput()+"/combineVCF."+String.valueOf(iter));
//            genotypedVariants.saveAsTextFile(options.getOutput()+"/genotypeVCF."+String.valueOf(iter));

            iter++;
            region.setStart(region.getEnd()+1);
            region.setEnd(region.getStart()+step);
            if(combineOutLocalFile.isDirectory()){
                String[] eles=combineOutLocalFile.list();
                for(String ele:eles){
                    File eleFile=new File(combineOutLocalFile.getAbsolutePath()+"/"+ele);
                    eleFile.delete();
                }
                combineOutLocalFile.delete();
            }

            //结束条件2：所有gvcf文件读到EOF
        }

        sc.stop();
    }

}
