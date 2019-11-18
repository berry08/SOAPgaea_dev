package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVCFOutputFormat;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.jointcalling.JointCallingEngine;
import org.bgi.flexlab.gaea.tools.jointcalling.util.MultipleVCFHeaderForJointCalling;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.VCFHeaderReader;
import org.seqdoop.hadoop_bam.util.WrapSeekable;

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class JointCallingReducer
		extends Reducer<WindowsBasedWritable, VariantContextWritable, NullWritable, VariantContextWritable> {

	private int windowSize;
	private HashMap<Integer, String> contigs = null;
	private JointCallingEngine engine = null;
	private VariantContextWritable outValue = new VariantContextWritable();

	private JointCallingOptions options = null;
	private GenomeLocationParser parser = null;
	private ReferenceShare genomeShare = null;
	private DbsnpShare dbsnpShare = null;
	private VCFLocalLoader loader = null;
	private VariantRegionFilter filter = null;
	private VCFHeader header = null;
	private MultipleVCFHeaderForJointCalling headers = new MultipleVCFHeaderForJointCalling();
	private ArrayList<ArrayList<String> > multiMapSampleNames=new ArrayList<ArrayList<String> >();
	SimpleDateFormat formatter = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss:SSS");
	@Override
	protected void setup(Context context) throws IOException {
		Configuration conf = context.getConfiguration();
		contigs = new HashMap<>();
		options = new JointCallingOptions();
		options.getOptionsFromHadoopConf(conf);
		Path path = new Path(conf.get(GaeaVCFOutputFormat.OUT_PATH_PROP));
		SeekableStream in = WrapSeekable.openPath(path.getFileSystem(conf), path);
//		SeekableStream in=new SeekableFileStream(new File(options.getTmpOut()+"/virtual.vcf"));
		System.out.println("before readHeader");
		header = VCFHeaderReader.readHeaderFrom(in);
		System.out.println("after readHeader");
		in.close();
		if(header == null)
			throw new RuntimeException("header is null !!!");
		
		for (VCFContigHeaderLine line : header.getContigLines()) {
			contigs.put(line.getContigIndex(), line.getID());
		}
		windowSize = options.getWindowsSize();
		parser = new GenomeLocationParser(header.getSequenceDictionary());
		headers.readHeaders(conf);
		String sampleStr = conf.get(JointCalling.INPUT_ORDER);;
		String[] allSample=sampleStr.split(",");
		int mapperLine=allSample.length/options.getMapperNumber()<2?2:allSample.length/options.getMapperNumber();
		ArrayList<String> ele=new ArrayList<String>();
		
		for(int i=0;i<allSample.length;i++) {
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
		System.out.println("before engine init");
		if(sampleStr != null) {
			if(options.getMapperMode() || options.getMapperMode3() || options.getMapperMode4()) {
				//JointCallingPrepareOptions options, GenomeLocationParser parser, VCFHeader vcfheader,
				//MultipleVCFHeaderForJointCalling multiHeaders,String[] sampleArray,ArrayList<ArrayList<String> >multiMapSamples,Configuration conf
				engine = new JointCallingEngine(options, parser,header,headers,allSample,multiMapSampleNames,conf);
			}else {
				engine = new JointCallingEngine(options, parser,headers);
			}
		}else {
			engine = new JointCallingEngine(options, parser,headers);
		}
		System.out.println("after engine init");
		genomeShare = new ReferenceShare();
		genomeShare.loadChromosomeList(options.getReference());
		dbsnpShare = new DbsnpShare(options.getDBSnp(), options.getReference());
		dbsnpShare.loadChromosomeList(options.getDBSnp() + VcfIndex.INDEX_SUFFIX);
		loader = new VCFLocalLoader(options.getDBSnp());
		filter = new VariantRegionFilter();
		System.out.println("reduce setup done");
	}

	@Override
	public void reduce(WindowsBasedWritable key, Iterable<VariantContextWritable> values, Context context)
			throws IOException, InterruptedException {
		System.out.println(formatter.format(new Date())+"\treduce start");
		int winNum = key.getWindowsNumber();
		int start = winNum * windowSize;
		if(start == 0)
			start = 1;
		String chr = contigs.get(key.getChromosomeIndex());
		int chrInx=key.getChromosomeIndex();
		int contigLength = header.getSequenceDictionary().getSequence(chr).getSequenceLength();
		int end = Math.min(contigLength, start + windowSize - 1);
		
		long startPosition = dbsnpShare.getStartPosition(chr, winNum, options.getWindowsSize());
		ArrayList<VariantContext> dbsnps = null;
		if(startPosition >= 0)
			dbsnps = filter.loadFilter(loader, chr, startPosition, end);
		engine.init(dbsnps);
		Set<Integer> bps=new TreeSet();
		String bpFile=options.getTmpOut()+"/AllBPs";
		BufferedReader bp_reader=new BufferedReader(new FileReader(bpFile));
		String winLine=null;
		while((winLine=bp_reader.readLine())!=null) {
			String[] bpeles=winLine.split("\t");
			Integer bpcontig=Integer.valueOf(bpeles[0]);
			Integer bpstart=Integer.valueOf(bpeles[1]);
			Integer bpend=Integer.valueOf(bpeles[2]);
			if(bpcontig==chrInx) {
				if(bpstart<=end && bpend>=start) {
					Integer realStart=bpstart>start?bpstart:start;
					Integer realEnd=bpend<end?bpend:end;
					for(int i=realStart;i<=realEnd;i++) {
						bps.add(i);
					}
				}else if(bpend<start) {
					continue;
				}else {
					break;
				}
			}else if(bpcontig<chrInx) {
				continue;
			}else {
				break;
			}
		}
		System.out.println(formatter.format(new Date())+"\tcurrent reduce key:\t"+chr+"\t"+start+"\t"+end);
//		for (int iter = start; iter <= end; iter++){
		for(Integer iter:bps) {
			VariantContext variantContext = engine.variantCalling(values.iterator(),
					parser.createGenomeLocation(chr, iter), genomeShare.getChromosomeInfo(chr));
			if (variantContext == null)
				continue;
			System.out.println(formatter.format(new Date())+"\tafter variantCalling\t"+variantContext.getChr()+"\t"+variantContext.getStart()+"\t"+variantContext.getEnd()+"\t"+variantContext.getReference()+"\t"+variantContext.getAlternateAlleles());
			CommonInfo info = variantContext.getCommonInfo();
			HashMap<String, Object> maps = new HashMap<>();
			maps.putAll(info.getAttributes());
			maps.remove("SM");
			info.setAttributes(maps);
			outValue.set(variantContext, header);
			context.write(NullWritable.get(), outValue);
		}
		//breakPoints.clear();
		//bp_reader.close();
	}
}
