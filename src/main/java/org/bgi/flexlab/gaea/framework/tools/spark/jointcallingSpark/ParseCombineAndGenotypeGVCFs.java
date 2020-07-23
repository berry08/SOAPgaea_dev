
package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.util.VariantContextHadoopWriter;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.annotator.util.Tuple;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.tools.mapreduce.jointcalling.JointCallingOptions;
import org.seqdoop.hadoop_bam.VariantContextCodec;
import org.seqdoop.hadoop_bam.VariantContextWithHeader;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import scala.Tuple2;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

public class ParseCombineAndGenotypeGVCFs implements Function2<Integer,Iterator<String>,Iterator<String>> {
    private HashMap<Integer, String> contigs = null;
    private JointCallingEngine engine = null;
    private HashMap<String,Integer> chrIndexs = new HashMap<>();
    private GenomeLocationParser parser = null;
    private ReferenceShare genomeShare = null;
    private DbsnpShare dbsnpShare = null;
    private VCFLocalLoader loader = null;
    private VariantRegionFilter filter = null;
    private VCFHeader header = null;
    private VCFHeader header2=null;
    private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
    private ArrayList<ArrayList<String>> multiMapSampleNames=new ArrayList<ArrayList<String> >();
    private VCFEncoder vcfEncoder=null;
    public BufferedReader bp_reader=null;
    public String winLine=null;
    public GenomeLongRegion region=null;
    SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
    HashMap<String, String> confMap=new HashMap<>();
    String bpFile;
    String[] args;
    DriverBC dBC;
    public  ArrayList<GenomeLocation> regions=new ArrayList<>();
    public int cycleIter;
    public ArrayList<Long> bpPartition=new ArrayList<>();
    public ParseCombineAndGenotypeGVCFs(GenomeLongRegion region, ArrayList<GenomeLocation> regions,String[] args, String outputBP, HashMap<String, String> confMap,
                            Broadcast<DriverBC> dBC,int cycleIter,ArrayList<Long> bpPartition) {
        this.args=args;
        this.dBC=dBC.value();
        this.cycleIter=cycleIter;
        this.confMap=confMap;
        this.region=region;
        this.regions=regions;
        bpFile=outputBP;
        for(Long ele:bpPartition){
            this.bpPartition.add(ele);
        }
    }
    public class VcComp implements Comparator<VariantContext>,Serializable{

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    Comparator<VariantContext> comparator4 = new VcComp();
    String HEADER_DEFAULT_PATH = "vcfheader";
    String MERGER_HEADER_INFO = "vcfheaderinfo";
    @Override public Iterator<String> call(Integer index,Iterator<String> stringIterator) throws IOException {
        Logger log= LoggerFactory.getLogger(ParseCombineAndGenotypeGVCFs.class);
        ArrayList<String> outValue=new ArrayList<>();
        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
        ArrayList<String> mapGvcfList=new ArrayList<String>();

        while(stringIterator.hasNext()) {
            String gvcfPath = stringIterator.next();
            mapGvcfList.add(gvcfPath);
        }
        contigs = new HashMap<>();
//        SeekableStream in=new SeekableFileStream(new File(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP)));
//		SeekableStream in=new SeekableFileStream(new File(options.getOutDir()()+"/virtual.vcf"));
        System.out.println(formatter.format(new Date())+"\tbefore readHeader");
        Path path = new Path(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
        SeekableStream in = new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        header = VCFHeaderReader.readHeaderFrom(in);
        System.out.println(formatter.format(new Date())+"\tafter readHeader");
        in.close();
        if(header == null)
            throw new RuntimeException("header is null !!!");

        for (VCFContigHeaderLine line : header.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
            chrIndexs.put(line.getID(),line.getContigIndex());
        }
        parser = new GenomeLocationParser(header.getSequenceDictionary());
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, dBC.options.getOutDir() + "/vcfheader");
        hadoop_conf.set(MERGER_HEADER_INFO, dBC.options.getOutDir()+"/"+MERGER_HEADER_INFO);
//        headers.readHeaders(hadoop_conf);
        SeekableStream in2=new SeekableFileStream(new File(dBC.outputDir+"/vcfheader"));
        headers.setMergeHeader(VCFHeaderReader.readHeaderFrom(in2));
        headers.setCurrentIndex(dBC.sampleIndex.size());
        HashMap<Integer, Set<String>> tmpNameHeaders=new HashMap<Integer, Set<String>>();
        for(Map.Entry<String,Integer> kv:dBC.sampleIndex.entrySet()){
            TreeSet<String> samples=new TreeSet<>();
            samples.add(kv.getKey());
            tmpNameHeaders.put(kv.getValue(),samples);
        }
        headers.setNameHeaders(tmpNameHeaders);
        String sampleStr = confMap.get(dBC.INPUT_ORDER);
        String[] allSample=sampleStr.split(",");
//        int mapperLine=allSample.length/options.getMapperNumber()<2?2:allSample.length/options.getMapperNumber();
        int mapperLine=10;
        ArrayList<String> ele=new ArrayList<String>();





        System.out.println(formatter.format(new Date())+"\tbefore engine init");
        dBC.options.setS(dBC.options.STANDARD_CONFIDENCE_FOR_CALLING);
        dBC.options.sets(dBC.options.STANDARD_CONFIDENCE_FOR_EMITTING);
        engine = new JointCallingEngine(dBC.options, parser,header,headers,allSample,dBC.multiMapSampleNames,hadoop_conf);
        System.out.println(formatter.format(new Date())+"\tafter engine init");
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(dBC.options.getReference());
        dbsnpShare = new DbsnpShare(dBC.options.getDBSnp(), dBC.options.getReference());
        dbsnpShare.loadChromosomeList(dBC.options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        loader = new VCFLocalLoader(dBC.options.getDBSnp());
        filter = new VariantRegionFilter();
        header2=engine.getVCFHeader();
        if(header2 == null)
            throw new RuntimeException("header is null !!!");
        Path pBPpath=new Path(bpFile);
        BufferedReader bp_reader=new BufferedReader(new FileReader(bpFile));
        winLine=bp_reader.readLine();
        System.out.println("reduce setup done");

        System.out.println(formatter.format(new Date())+"\treduce start");
        ArrayList<BufferedReader> samplesReader=new ArrayList<>();
        ArrayList<VCFCodec> samplesCodec=new ArrayList<>();
        ArrayList<Integer> samplesTag=new ArrayList<>();
        for(ArrayList<String> samplesInMap:dBC.multiMapSampleNames){
            //打开文件句柄
            Integer mapSMtagInt=0;
            int samplesInMapIndex=0;
            for(String sampleName:samplesInMap){
                if(samplesInMapIndex==0) {
                    mapSMtagInt = dBC.sampleIndex.get(sampleName);
                }
                samplesInMapIndex++;
            }
            mapSMtagInt+=headers.getHeaderSize();
            samplesTag.add(mapSMtagInt);
//            String targetFileName=dBC.options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt+"."+index;
            String targetFileName=dBC.options.getOutDir()+"/combine."+cycleIter+"/combineFile."+mapSMtagInt+"."+index;
            System.out.println("processed file:\t"+targetFileName);
//            Path filePath=new Path(targetFileName);
//            FileSystem fileFs=filePath.getFileSystem(hadoop_conf);
//            if(!fileFs.exists(filePath)){
//                log.error("code error1\t"+filePath.getName()+" not exists");
//                System.exit(1);
//            }

//            BufferedReader reader=new BufferedReader(new InputStreamReader(fileFs.open(filePath)));
            BufferedReader reader=new BufferedReader(new FileReader(targetFileName));
            VCFHeader curSamplesMergedHeader=new VCFHeader(dBC.virtualHeader.getMetaDataInInputOrder(),samplesInMap);
            VCFCodec codec=new VCFCodec();
            codec.setVCFHeader(curSamplesMergedHeader,dBC.version);
            samplesCodec.add(codec);
            samplesReader.add(reader);
        }
        if(samplesReader.size()!=dBC.multiMapSampleNames.size() || samplesCodec.size()!=dBC.multiMapSampleNames.size()){
            log.error("code error2\t"+samplesReader.size()+"\t"+dBC.multiMapSampleNames.size());
            System.exit(1);
        }
        int vcNum=0;
        String outFile=dBC.options.getOutDir()+"/genotype."+cycleIter+"/"+region.getStart()+"_"+region.getEnd()+"."+index;
        Path outPath=new Path(outFile);
        FileSystem outFs=outPath.getFileSystem(hadoop_conf);
        if(outFs.exists(outPath)){
            outFs.delete(outPath,true);
        }
        BufferedWriter out=new BufferedWriter(new FileWriter(outFile));
        long rangeInEachOutFile=(region.getEnd()-region.getStart())/dBC.options.getReducerNumber();
        log.info("range in each out file:\t"+rangeInEachOutFile);
        log.info("current index:\t"+index);
        int minRange=100000;
        if(rangeInEachOutFile<minRange){
            rangeInEachOutFile=minRange;
        }
        boolean[] readerDone=new boolean[dBC.multiMapSampleNames.size()];
        for(int i=0;i<dBC.multiMapSampleNames.size();i++){
            readerDone[i]=false;
        }
        vcfEncoder = new VCFEncoder(header, true, true);
        for(GenomeLocation curRegion:regions) {
            int start = curRegion.getStart();
            String chr = curRegion.getContig();
            int chrInx = dBC.chrIndex.get(chr);
            int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
            int end = curRegion.getEnd();
            long regionStart=dBC.accumulateLength.get(chrInx)+start;
            long regionEnd=dBC.accumulateLength.get(chrInx)+end;
            if(regionStart>bpPartition.get(index)){
                break;
            }else{
                if(index>0){
                    if(regionEnd<bpPartition.get(index-1)){
                        continue;
                    }
                }
            }
            ArrayList<VariantContext> dbsnps = null;
            long startPosition = dbsnpShare.getStartPosition(chr, start / dBC.options.getWindowsSize(), dBC.options.getWindowsSize());
            if (startPosition >= 0) dbsnps = filter.loadFilter(loader, chr, startPosition, end);
            engine.init(dbsnps);



            //为了避免每次从头读取winBp，保留bp_reader，使其只读一遍
            System.out.println(formatter.format(new Date()) + "\tbp window before get bps:\t" + winLine);
            LinkedList<VariantContext> regionVcs=new LinkedList<>();
//            long smallWinStartLong=dBC.accumulateLength.get(chrInx)+start+rangeInEachOutFile*index;
            long smallWinStartLong=regionStart;
            long smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
            int smallWinStart=start;
            int smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
            HashMap<Integer,ArrayList<VariantContext>> lastVCs=new HashMap<>();
            ArrayList<BufferedWriter> outWriterList=new ArrayList<>();
            int logIter=0;
            while(smallWinStartLong<=regionEnd){
                //一个个小窗口开始处理
                //获取小窗口内的所有变异
                if(smallWinStartLong>bpPartition.get(index)){
                    break;
                }
                if(index>0){
                    if(smallWinEndLong<bpPartition.get(index-1)) {
                        smallWinStartLong=smallWinEndLong+1;
                        smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
                        smallWinStart=smallWinEnd+1;
                        smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
                        continue;
                    }
                }
//                int part1= (int) ((smallWinStartLong-region.getStart())/rangeInEachOutFile);
//                int part2= (int) ((smallWinEndLong-region.getStart())/rangeInEachOutFile);
//                if(part2<index){
//                    smallWinStartLong=smallWinEndLong+1;
//                    smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
//                    smallWinStart=smallWinEnd+1;
//                    smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
//                    continue;
//                }else if(part1>index){
//                    break;
//                }

                Set<Integer> bps=new TreeSet<>();
                for(int i=0;i<dBC.multiMapSampleNames.size();i++){
                    if(readerDone[i]){
                        continue;
                    }
                    boolean remove=false;
                    if(lastVCs.containsKey(i)){
                        for(VariantContext lastVC:lastVCs.get(i)){
                            long curPosStartLong=dBC.accumulateLength.get(dBC.chrIndex.get(lastVC.getContig()))+lastVC.getStart();
                            long curPosEndLong=dBC.accumulateLength.get(dBC.chrIndex.get(lastVC.getContig()))+lastVC.getEnd();
                            if(curPosStartLong<=smallWinEndLong && curPosEndLong>=smallWinStartLong){
                                remove=true;
                                if(lastVC.getNAlleles()>2) {
                                    for (int pos = lastVC.getStart(); pos <= lastVC.getEnd(); pos++) {
                                        if(pos>=smallWinStart && pos<=smallWinEnd)
                                            bps.add(pos);
                                    }
                                }
                                regionVcs.add(lastVC);
                            }else if(curPosEndLong<smallWinStartLong){
                                remove=true;
                            }

                        }
                        if(remove){
                            lastVCs.remove(i);
                        }
                    }
                    String line = null;
                    while(true){
                        try {
                            line=samplesReader.get(i).readLine();
                            if(line==null){
                                readerDone[i]=true;
                                break;
                            }
                        } catch (IOException e) {
                            System.out.println("IO error!\tcurrent sample:\t"+samplesTag+"\t"+cycleIter);
                            System.out.println("current window:\t"+smallWinStartLong);
                            if(lastVCs.containsKey(i)) {
                                System.out.println("last VC:");
                                for(VariantContext lastVC:lastVCs.get(i)){
                                    System.out.println(lastVC);
                                }
                            }
                            e.printStackTrace();
                        }

                        VariantContext combineVC=samplesCodec.get(i).decode(line);
                        long curPosStartLong=dBC.accumulateLength.get(dBC.chrIndex.get(combineVC.getContig()))+combineVC.getStart();
                        long curPosEndLong=dBC.accumulateLength.get(dBC.chrIndex.get(combineVC.getContig()))+combineVC.getEnd();
                        if(curPosStartLong<=smallWinEndLong && curPosEndLong>=smallWinStartLong){
                            if(combineVC.getNAlleles()>2) {
                                for (int pos = combineVC.getStart(); pos <= combineVC.getEnd(); pos++) {
                                    if(pos>=smallWinStart && pos<=smallWinEnd)
                                        bps.add(pos);
                                }
                            }
                            regionVcs.add(combineVC);
                            if(curPosEndLong>smallWinEndLong){
                                if(lastVCs.containsKey(i)){
                                    lastVCs.get(i).add(combineVC);
                                }else{
                                    ArrayList<VariantContext> tmpVCs=new ArrayList<>();
                                    tmpVCs.add(combineVC);
                                    lastVCs.put(i,tmpVCs);
                                }
                                break;
                            }
                        }else if(curPosStartLong>smallWinEndLong){
                            if(lastVCs.containsKey(i)){
                                lastVCs.get(i).add(combineVC);
                            }else{
                                ArrayList<VariantContext> tmpVCs=new ArrayList<>();
                                tmpVCs.add(combineVC);
                                lastVCs.put(i,tmpVCs);
                            }
                            break;
                        }else{
                            continue;
                        }
                    }
                }
                Collections.sort(regionVcs, comparator4);
                if(logIter%5==0) {
                    System.out.println(formatter.format(new Date()) + "\tcurrent reduce key:\t" + chr + "\t" + chrInx + "\t" + smallWinStartLong + "\t" + smallWinEndLong);
                    System.out.println(formatter.format(new Date()) + "\tbps size in region:\t" + bps.size());
                }
                logIter++;
                for(int iter:bps){
                    VariantContext variantContext = engine.variantCallingForSpark(regionVcs.iterator(),
                            parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
                    if (variantContext == null){
                        continue;
                    }
                    CommonInfo info = variantContext.getCommonInfo();
                    HashMap<String, Object> maps = new HashMap<>();
                    maps.putAll(info.getAttributes());
                    maps.remove("SM");
                    info.setAttributes(maps);
                    String value=vcfEncoder.encode(variantContext);
                    out.write(value+"\n");
                    vcNum++;
                    while(regionVcs.size()>0){
                        VariantContext firstVc=regionVcs.getFirst();
                        if(firstVc.getEnd()<=iter){
                            regionVcs.removeFirst();
                        }else{
                            break;
                        }
                    }
                }
                smallWinStartLong=smallWinEndLong+1;
                smallWinEndLong=smallWinStartLong+dBC.options.getWindowsSize()-1;
                smallWinStart=smallWinEnd+1;
                smallWinEnd=smallWinStart+dBC.options.getWindowsSize()-1;
                regionVcs.clear();
            }
        }
        out.close();
        ArrayList<String> totalNum=new ArrayList<>();
        totalNum.add("done"+index);
        for(int i=0;i<samplesReader.size();i++){
            samplesReader.get(i).close();
        }
        return totalNum.iterator();
    }
}
