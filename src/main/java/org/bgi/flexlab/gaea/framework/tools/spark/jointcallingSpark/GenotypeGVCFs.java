package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.VoidFunction;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.VCFHdfsWriter;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;


public class GenotypeGVCFs implements VoidFunction<Iterator<String>> {
    public JointCallingSparkOptions options;
    public GenotypeGVCFs(JointCallingSparkOptions options){
        this.options=options;
    }
    public class VcComp implements Comparator<VariantContext>, Serializable {

        @Override public int compare(VariantContext o1, VariantContext o2) {
            return o1.getStart()-o2.getStart();
        }
    }
    Comparator<VariantContext> comparator4 = new GenotypeGVCFs.VcComp();
    @Override public void call(Iterator<String> stringIterator) throws Exception {
//        Runtime runtime = Runtime.getRuntime();
////        Process pro=runtime.exec("ulimit -n 7654321");
////        int status = pro.waitFor();
////        if (status != 0)
////        {
////            System.out.println("Failed to call shell's command ");
////            System.exit(1);
////        }
////        pro=runtime.exec("ulimit -a");
        Configuration conf=new Configuration();
        String MERGER_HEADER_INFO = "vcfheaderinfo";
        conf.set(GaeaVCFOutputFormat.OUT_PATH_PROP, options.getVCFHeaderOutput() + "/vcfheader");
        conf.set(MERGER_HEADER_INFO, options.getVCFHeaderOutput()+"/"+MERGER_HEADER_INFO);
        List<String> gvcfList=options.getInputStringList();
        BufferedReader[] samplesReader=new BufferedReader[gvcfList.size()];
        Logger log= LoggerFactory.getLogger(GenotypeGVCFs.class);
        HashMap<String,Integer> sampleIndex=new HashMap<>();
        String sampleIndexFile=options.getOutDir()+"/vcfheaderinfo";
        Path sampleIndexFilePath=new Path(sampleIndexFile);
        BufferedReader vcfheaderinfoReader= new BufferedReader(new InputStreamReader(sampleIndexFilePath.getFileSystem(conf).open(sampleIndexFilePath)));
        String tmpLine;
        while((tmpLine=vcfheaderinfoReader.readLine())!=null){
            String[] eles = tmpLine.split("\t");
            if (eles.length != 3) {
                log.error("vcfheaderinfo file format error");
            }
            String name;
            if (eles[2].endsWith(",")) {
                name = eles[2].substring(0, eles[2].length() - 1);
            } else {
                name = eles[2];
            }
            sampleIndex.put(name,Integer.valueOf(eles[1]));
//            pathSample.put(eles[0], name);
        }
        int i=0;
        log.info("before open files");
        int bufsize=1024;
        for(String gvcfSamplePath:gvcfList){
            if(gvcfSamplePath.startsWith("file://")) {
                gvcfSamplePath=gvcfSamplePath.substring(7);
            }
            BufferedReader reader=null;
            if(gvcfSamplePath.endsWith(".gz")) {
                if(gvcfSamplePath.startsWith("/user")) {
                    Path samplePath=new Path(gvcfSamplePath);
                    FileSystem gvcf_fs=samplePath.getFileSystem(conf);
                    reader=new BufferedReader(new InputStreamReader(new GZIPInputStream(gvcf_fs.open(new Path(gvcfSamplePath)))),bufsize);
                }else {
                    reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gvcfSamplePath))),bufsize);
                }
            }else {
                reader = new BufferedReader(new InputStreamReader(new FileInputStream(gvcfSamplePath)),bufsize);
            }
            samplesReader[i]=reader;
            i++;
        }
        log.info("after open files");
        log.info("open files number:\t"+samplesReader.length);

        Path mergedHeaderPath=new Path(options.getOutDir()+"/vcfheader");
        SeekableStream in = WrapSeekable.openPath(mergedHeaderPath.getFileSystem(conf), mergedHeaderPath);
        VCFHeader header = VCFHeaderReader.readHeaderFrom(in);
        //生成codec
        Path virualHeaderPath=new Path(options.getOutDir()+"/virtual.vcf");
        SeekableStream in2 = WrapSeekable.openPath(virualHeaderPath.getFileSystem(conf), virualHeaderPath);
        VCFHeader virtualHeader = VCFHeaderReader.readHeaderFrom(in2);
        int iter=0;
        Set<VCFHeaderLine> gvcfHeaderMetaInfo=virtualHeader.getMetaDataInInputOrder();
        VCFHeaderVersion version=null;
        HashMap<String,Integer> chrIndex=new HashMap<>();
        HashMap<Integer,String> contigs=new HashMap<>();
        for (VCFContigHeaderLine line : virtualHeader.getContigLines()) {
            contigs.put(line.getContigIndex(), line.getID());
            chrIndex.put(line.getID(),line.getContigIndex());
        }
        for (final VCFHeaderLine line : gvcfHeaderMetaInfo) {
            if (VCFHeaderVersion.isFormatString(line.getKey())) {
                version = VCFHeaderVersion.toHeaderVersion(line.getValue());
                break;
            }
        }
        VCFHdfsWriter vcfHdfsWriter=null;

        VCFCodec[] samplesCodec=new VCFCodec[samplesReader.length];

        MultipleVCFHeaderForJointCalling headers=new MultipleVCFHeaderForJointCalling();
        headers.readHeaders(conf);
        GenomeLocationParser parser=new GenomeLocationParser(virtualHeader.getSequenceDictionary());
        options.setS(options.STANDARD_CONFIDENCE_FOR_CALLING);
        options.sets(options.STANDARD_CONFIDENCE_FOR_EMITTING);
        JointCallingEngine engine = new JointCallingEngine(options, parser,headers);
        ReferenceShare genomeShare = new ReferenceShare();
        genomeShare.loadChromosomeList(options.getReference());
        VCFLocalLoader loader = new VCFLocalLoader(options.getDBSnp());
        VariantRegionFilter filter = new VariantRegionFilter();

        String[] sampleNames=new String[samplesReader.length];
        for(int ii=0;ii<samplesReader.length;ii++){
            String line;
            while(true){
                line=samplesReader[ii].readLine();
                if(line.startsWith("##")){
                    continue;
                }else if(line.startsWith("#")){
                    String[] eles=line.split("\t");
                    sampleNames[ii]=eles[eles.length-1];
                    VCFCodec codec=new VCFCodec();
                    ArrayList<String> curSamples=new ArrayList<>();
                    curSamples.add(eles[eles.length-1]);
                    virtualHeader=new VCFHeader(gvcfHeaderMetaInfo,curSamples);
                    codec.setVCFHeader(virtualHeader,version);
                    samplesCodec[ii]=codec;
                    break;
                }
            }
        }
        VCFCodec codec=new VCFCodec();
        codec.setVCFHeader(virtualHeader,version);
//        log.info("before seek");
//        long skipSize=1024*1024*100;
//        for(int ii=0;ii<samplesReader.length;ii++){
//            String line2;
//            while((line2=samplesReader[ii].readLine())!=null){
//                if((line2=samplesReader[ii].readLine())!=null){
//                    log.info("sample "+ii+"\t"+line2);
//                    samplesReader[ii].skip(skipSize);
//                }
//            }
//            if(ii%1000==0){
//                log.info("sample "+i+" seek done");
//            }
//        }
//
//        log.info("after seek");
//        return;
        VariantContext[] curSamplesVC=new VariantContext[samplesReader.length];
        for(int ii=0;ii<samplesReader.length;ii++){
            curSamplesVC[ii]=null;
        }
        ArrayList<VariantContext> dbsnps=new ArrayList<>();
        VCFEncoder vcfEncoder = new VCFEncoder(header, true, true);
        BufferedWriter out=null;
        log.info("start circulation");
        long[] accumulateLength=new long[chrIndex.size()];
        long totalLength=0;
        for(int ii=0;ii<chrIndex.size();ii++){
            long chrLength=header.getSequenceDictionary().getSequence(ii).getSequenceLength();
            totalLength+=chrLength;
            accumulateLength[ii]=totalLength-chrLength;
        }
        while(stringIterator.hasNext()){
            String region=stringIterator.next();
            log.info("process region:\t"+region);
            String[] eles=region.split("\t");
            if(eles.length!=3){
                log.error("region parse error,"+region);
                System.exit(1);
            }
            String chr=eles[0];
            Integer chrInt=chrIndex.get(chr);
            Integer start=Integer.parseInt(eles[1]);
            Integer end=Integer.parseInt(eles[2]);


            if(iter==0){
                //创建输出文件
                vcfHdfsWriter = new VCFHdfsWriter(options.getOutput()+"/"+chr+"_"+start, false, false, conf);
                Path outPath=new Path(options.getOutput()+"/"+chr+"_"+start);
                FileSystem fs=outPath.getFileSystem(conf);
                vcfHdfsWriter.writeHeader(header);//写头文件
                vcfHdfsWriter.close();
//                out=new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(fs.create(outPath))));
                out=new BufferedWriter(new OutputStreamWriter(fs.create(outPath)));
//                out.write(header.toString()+"\n");
                log.info("before seek");
                long startLong=accumulateLength[chrInt]+start;
                int indexPos= (int) (startLong/1000000);
                if(indexPos>0) {
                    Path indexPath = new Path(options.getOutDir() + "/index/" + indexPos);
                    FileSystem indexFs = indexPath.getFileSystem(conf);
                    if (!indexFs.exists(indexPath)) {
                        log.info("code error");
                        System.exit(1);
                    }
                    BufferedReader simpleIndexReader = new BufferedReader(new InputStreamReader(indexFs.open(indexPath)));
                    String line;
                    HashMap<Integer, Long> samplesOffset = new HashMap<>();
                    while ((line = simpleIndexReader.readLine()) != null) {
                        String[] eles2 = line.split("\t");
                        long pos = Long.parseLong(eles2[1]);
                        String samplename = eles2[0];
                        samplesOffset.put(sampleIndex.get(samplename), pos);
                    }
                    simpleIndexReader.close();
                    for (int ii = 0; ii < samplesReader.length; ii++) {
                        samplesReader[ii].skip(samplesOffset.get(ii));
                        log.info("after search. Target offset is " + samplesOffset.get(ii) + ". Sample " + ii + " done");
                    }
                }
                log.info("after seek");
            }


            //处理一段区域，variantCalling后写入vcf
            ArrayList<VariantContext> regionVcs=new ArrayList<>();
            DbsnpShare dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
            dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
            long startPosition = dbsnpShare.getStartPosition(chr, start/options.getWindowsSize(), options.getWindowsSize());
            if(startPosition >= 0) {
                dbsnps = filter.loadFilter(loader, chr, startPosition, end);
            }
            engine.init(dbsnps);
            Set<Integer> bps=new TreeSet<>();
            log.info("before query samples");

            for(int ii=0;ii<samplesReader.length;ii++){
                //处理上一次
                if(curSamplesVC[ii]!=null){
                    if(curSamplesVC[ii].getContig().equals(chr) && curSamplesVC[ii].getStart()<=end && curSamplesVC[ii].getEnd()>=start){
                        CommonInfo info = curSamplesVC[ii].getCommonInfo();
                        if(!info.hasAttribute("SM")) {
                            if(sampleIndex.containsKey(sampleNames[ii])) {
                                info.putAttribute("SM", sampleIndex.get(sampleNames[ii]));
                            }else{
                                log.error("no such sample index");
                            }
                        }
                        if(curSamplesVC[ii].getNAlleles()>2){
                            for(int s=curSamplesVC[ii].getStart();s<=curSamplesVC[ii].getEnd();s++){
                                bps.add(s);
                            }
                        }
                        regionVcs.add(curSamplesVC[ii]);
                    }
                    curSamplesVC[ii]=null;
                }
                ArrayList<String> curSampleName=new ArrayList<>();
                curSampleName.add(sampleNames[ii]);

                String line;
                int iter2=0;
                while((line=samplesReader[ii].readLine())!=null){
                    VariantContext v=samplesCodec[ii].decode(line);
                    iter2++;
//                    VariantContext v=codec.decode(line);
                    if(iter==0 && iter2<2){
                        log.info("region:\t"+region+"\tsample "+ii+"\tfirst variant:"+line);
                    }
                    if(v.getContig().equals(chr)){
                        if(v.getStart()<=end && v.getEnd()>=start){
                            CommonInfo info = v.getCommonInfo();
                            GenotypesContext gc=v.getGenotypes();
                            if(!info.hasAttribute("SM")) {
                                if(sampleIndex.containsKey(sampleNames[ii])) {
                                    info.putAttribute("SM", sampleIndex.get(sampleNames[ii]));
                                }else{
                                    log.error("no such sample index");
                                }
                            }
//                            FastGenotype tmpGc=new FastGenotype(sampleNames[ii],v.getAlleles(),false,gc.get(0).getGQ(),gc.get(0).getDP(),gc.get(0).getAD(),gc.get(0).getPL(),gc.get(0).getFilters(),gc.get(0).getExtendedAttributes());
                            regionVcs.add(v);
                            if(v.getNAlleles()>2){
                                for(int s=v.getStart();s<=v.getEnd();s++){
                                    bps.add(s);
                                }
                            }
                            if(v.getEnd()>end){
                                curSamplesVC[ii]=v;
                                break;
                            }
                        }else{
                            if(v.getStart()>end){
                                curSamplesVC[ii]=v;
                                break;
                            }else{
                                //可以快速跳转
                                continue;
                            }
                        }
                    }else{
                        if(chrIndex.get(v.getContig())>chrInt){
                            curSamplesVC[ii]=v;
                            break;
                        }else{
                            //可以快速跳转
                            continue;
                        }
                    }
                }
            }
            log.info("after query samples");
            log.info("before calling. VC number in this region:\t"+regionVcs.size()+"\tbp number:\t"+bps.size());
            Collections.sort(regionVcs,comparator4);
            for(Integer ii:bps) {
                VariantContext variantContext = engine.variantCallingForSpark(regionVcs.iterator(), parser.createGenomeLocation(chr, ii), genomeShare.getChromosomeInfo(chr));
                if (variantContext == null) {
                    continue;
                }
                CommonInfo info = variantContext.getCommonInfo();
                HashMap<String, Object> maps = new HashMap<>();
                maps.putAll(info.getAttributes());
                maps.remove("SM");
                info.setAttributes(maps);
                out.write(vcfEncoder.encode(variantContext)+"\n");
            }
            log.info("after calling");
//            out.flush();
            iter++;
        }
        out.close();
        log.info("end circulation");
    }

}
