package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.Function2;
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
import scala.Tuple2;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

public class GenotypeVariantsNoMerge implements Function2<Integer,Iterator<VariantContext>,Iterator<String>> {
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
    public GenomeLocation region=null;
    SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
    HashMap<String, String> confMap=new HashMap<>();
    String bpFile;
    String[] args;
    DriverBC dBC;
    public GenotypeVariantsNoMerge(GenomeLocation region, String[] args, String outputBP, HashMap<String, String> confMap,
                            Broadcast<DriverBC> dBC) {
        this.args=args;
        this.dBC=dBC.value();
        this.confMap=confMap;
        this.region=region;
        bpFile=outputBP;
    }
    String HEADER_DEFAULT_PATH = "vcfheader";
    String MERGER_HEADER_INFO = "vcfheaderinfo";
    @Override public Iterator<String> call(Integer v1,Iterator<VariantContext> tupleIterator) throws
            Exception {
        ArrayList<String> outValue=new ArrayList<>();
        JointCallingOptions options = new JointCallingOptions();
        options.parse(args);
        Configuration hadoop_conf=new Configuration();
        for(String k:confMap.keySet()){
            hadoop_conf.set(k,confMap.get(k));
        }
        contigs = new HashMap<>();
//        SeekableStream in=new SeekableFileStream(new File(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP)));
//		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
        System.out.println(formatter.format(new Date())+"\tbefore readHeader");
        Path path = new Path(confMap.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
        SeekableStream in = WrapSeekable.openPath(path.getFileSystem(hadoop_conf), path);
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
        hadoop_conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/vcfheader");
        hadoop_conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutput()+"/"+MERGER_HEADER_INFO);
        headers.readHeaders(hadoop_conf);

        String sampleStr = confMap.get(dBC.INPUT_ORDER);
        String[] allSample=sampleStr.split(",");
//        int mapperLine=allSample.length/options.getMapperNumber()<2?2:allSample.length/options.getMapperNumber();
        int mapperLine=10;
        ArrayList<String> ele=new ArrayList<String>();

        for(int i=0;i<allSample.length;i++) {
            String[] eles=allSample[i].split("/");
            Path path2=new Path(allSample[i]);
//            if(!pathSample.containsKey(eles[eles.length-1])) {
//                logger.error("no such path in vcfHeaderInfo");
//            }
//            String sampleName=pathSample.get(eles[eles.length-1]);
//            sampleNames.add(sampleName);
            ele.add(allSample[i]);
            if((i+1)%mapperLine==0) {
                ArrayList<String> tmp_ele=new ArrayList<String>();
                tmp_ele.addAll(ele);
                Collections.sort(tmp_ele);
                multiMapSampleNames.add(tmp_ele);
                ele.clear();
            }
        }
        if(ele.size()>0) {
            ArrayList<String> tmp_ele=new ArrayList<String>();
            tmp_ele.addAll(ele);
            Collections.sort(tmp_ele);
            multiMapSampleNames.add(tmp_ele);
        }
        System.out.println(formatter.format(new Date())+"\tbefore engine init");
//        engine = new JointCallingEngine(options, parser,header,headers,allSample,dBC.multiMapSampleNames,hadoop_conf);
        engine = new JointCallingEngine(options, parser,headers);
        System.out.println(formatter.format(new Date())+"\tafter engine init");
        genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());
        dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
        dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
        loader = new VCFLocalLoader(options.getDBSnp());
        filter = new VariantRegionFilter();
        header2=engine.getVCFHeader();
        if(header2 == null)
            throw new RuntimeException("header is null !!!");
        Path pBPpath=new Path(bpFile);
        BufferedReader bp_reader=new BufferedReader(new InputStreamReader(pBPpath.getFileSystem(hadoop_conf).open(pBPpath)));
        winLine=bp_reader.readLine();
        System.out.println("reduce setup done");

        System.out.println(formatter.format(new Date())+"\treduce start");
        int start = region.getStart();
        String chr = region.getContig();
        int chrInx=region.getContigIndex();
        int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
        int end = region.getEnd();
        int band=(region.getEnd()-region.getStart())/dBC.options.getReducerNumber();
        if(v1==0){
            start=region.getStart();
            end=start+band;
        }else{
            start=region.getStart()+band*v1;
            end=start+band;
        }


        ArrayList<VariantContext> dbsnps = null;
        long startPosition = dbsnpShare.getStartPosition(chr, start/options.getWindowsSize(), options.getWindowsSize());
        if(startPosition >= 0)
            dbsnps = filter.loadFilter(loader, chr, startPosition, end);
        engine.init(dbsnps);
        Set<Integer> bps=new TreeSet();

        vcfEncoder=new VCFEncoder(header, true, true);

        //为了避免每次从头读取winBp，保留bp_reader，使其只读一遍
//        System.out.println(formatter.format(new Date())+"\tbp window before get bps:\t"+winLine);
//        while(true) {
//            String[] bpeles=winLine.split(":");
//            Integer bpcontig=chrIndexs.get(bpeles[0]);
//            Integer bpstart=0;
//            Integer bpend=0;
//            if(bpeles[1].contains("-")){
//                String[] pairPos=bpeles[1].split("-");
//                bpstart=Integer.valueOf(pairPos[0]);
//                bpend=Integer.valueOf(pairPos[1]);
//            }else{
//                bpend=bpstart=Integer.valueOf(bpeles[1]);
//            }
//            if(bpcontig==chrInx) {
//                if(bpstart<=end && bpend>=start) {
//                    Integer realStart=bpstart>start?bpstart:start;
//                    Integer realEnd=bpend<end?bpend:end;
//                    for(int i=realStart;i<=realEnd;i++) {
//                        bps.add(i);
//                    }
//                    if(bpend<end) {
//                        if ((winLine = bp_reader.readLine()) == null) {
//                            break;
//                        }
//                    }else{
//                        break;
//                    }
//                }else if(bpend<start) {
//                    if((winLine=bp_reader.readLine())==null){
//                        break;
//                    }
//                }else if(bpstart>end){
//                    break;
//                }
//            }else if(bpcontig<chrInx) {
//                if((winLine=bp_reader.readLine())==null){
//                    break;
//                }
//            }else {
//                break;
//            }
//        }
//        System.out.println(formatter.format(new Date())+"\tcurrent reduce key:\t"+chr+"\t"+chrInx+"\t"+start+"\t"+end);
//        System.out.println(formatter.format(new Date())+"\tbp window after get bps:\t"+winLine);
//        System.out.println(formatter.format(new Date())+"\tbps size in region:\t"+bps.size());
//		for (int iter = start; iter <= end; iter++){
        int iterTmp=0;
//        while(tupleIterator.hasNext()){
//            System.out.println(tupleIterator.next());
//        }

//        for(Integer iter:bps) {
        for(Integer iter=start;iter<=end;iter++){
            VariantContext variantContext = engine.variantCallingForSpark2(tupleIterator,
                    parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
            if (variantContext == null){
                continue;
            }
//			iterTmp++;
//			if(iterTmp<3){
//				System.out.println(variantContext);
//				System.out.println(variantContext.getAttribute("BaseQRankSum"));
//			}
            CommonInfo info = variantContext.getCommonInfo();
            HashMap<String, Object> maps = new HashMap<>();
            maps.putAll(info.getAttributes());
            maps.remove("SM");
            info.setAttributes(maps);
//            VariantContextWithHeader vcWithHeader=new VariantContextWithHeader(variantContext,header);
            String value=vcfEncoder.encode(variantContext);
            outValue.add(value);
//            outValue.set(variantContext, header2);
//            context.write(NullWritable.get(), outValue);
        }
        return outValue.iterator();
    }
}
