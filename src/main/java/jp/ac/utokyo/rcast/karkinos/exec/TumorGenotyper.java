/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.rcast.karkinos.exec;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.karkinos.ploidy.MatchMatrixBean;
import jp.ac.utokyo.karkinos.ploidy.PloidyResolve;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.CNVUtils;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.CheckPossibleHDAmp;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.LoadSave;
import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.SaveBean;
import jp.ac.utokyo.rcast.karkinos.filter.DefinedSites;
import jp.ac.utokyo.rcast.karkinos.filter.FilterAnnotation;
import jp.ac.utokyo.rcast.karkinos.filter.SoftClipExtention;
import jp.ac.utokyo.rcast.karkinos.filter.TerminalMismatch;
import jp.ac.utokyo.rcast.karkinos.graph.output.CNVVcf;
import jp.ac.utokyo.rcast.karkinos.graph.output.FileOutPut;
import jp.ac.utokyo.rcast.karkinos.graph.output.PdfReport;
import jp.ac.utokyo.rcast.karkinos.graph.output.TextSummary;
import jp.ac.utokyo.rcast.karkinos.hmm.CountCNV;
import jp.ac.utokyo.rcast.karkinos.hmm.HMMCNVAnalysisFromEM;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.utils.CorrelVaridate;
import jp.ac.utokyo.rcast.karkinos.utils.GeneEachCNV;
import jp.ac.utokyo.rcast.karkinos.utils.Interval;
import jp.ac.utokyo.rcast.karkinos.utils.ListUtils;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.SamUtils;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinos.wavelet.EMMethod;
import jp.ac.utokyo.rcast.karkinos.wavelet.GCParcentAdjust;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletDenoize;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class TumorGenotyper extends ReadWriteBase {
	public static void main(String[] arg) {
		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOptionListForKarkinos();
		Options opts = new Options();
		for (Option opt : optionList) {
			opts.addOption(opt);
		}

		CommandLine cl = null;
		try {
			cl = parcer.parse(opts, arg);
		} catch (ParseException e1) {
			System.out.println(e1.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.setOptionComparator(new OptionComparator(optionList));
			help.printHelp("karkinos.jar analysis", opts, true);
			return;
		}

		TumorGenotyper tg = new TumorGenotyper();
		try {
			tg.exec(cl);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static List<Option> getOptionListForKarkinos() {
		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("n", "normalBam", true, "normal bam file",
				true));
		optionlist
				.add(getOption("t", "tumorBam", true, "tumor bam file", true));
		optionlist.add(getOption("r", "reference", true,
				"2 bit genome reference file", true));
		optionlist.add(getOption("snp", "dbSNP", true,
				"dbSNP list from annover sites,(bin,chr,start,end)", true));
		optionlist.add(getOption("ct", "captureTarget", true,
				"Capture target regions(bed format)", false));

		optionlist
				.add(getOption("o", "outdir", true, "output directory", true));
		optionlist.add(getOption("id", "uniqueid", true,
				"unique id for this sample", true));
		optionlist.add(getOption("prop", "property", true,
				"path to property file( otherwise default val)", false));
		optionlist.add(getOption("mp", "mappability", true,
				"optional,mappability from ucsc (bw, big wig format)", false));
		optionlist
				.add(getOption(
						"g1000",
						"1000genome",
						true,
						"optional,1000 genome list from annover sites,(chr,pos,ref,alt,freq,id)",
						false));

		optionlist.add(getOption("cosmic", "cosmicSNV", true,
				"cosmic snv vcf format", false));

		optionlist.add(getOption("g1000freq", "1000genomefreq", true,
				"optional,1000 genome frequency threshold to use", false));

		optionlist.add(getOption("chr", "chrom", true, "chromosome to analyze",
				false));

		optionlist.add(getOption("startend", "startend", true,
				"start-end position", false));

		optionlist.add(getOption("rs", "readsStats", true,
				"optional,reads stats files(normal,tumor)", false));

		optionlist.add(getOption("rg", "refFlatGenes", true,
				"optional,gene reference for depth stats", false));

		optionlist.add(getOption("exonSNP", "exonSNP", true,
				"additional Exon SNP", false));

		optionlist.add(getOption("nopdf", "nopdf", false,
				"no graphic summary pdf output", false));


		optionlist.add(getOption("sites", "pileupsites", true,
				"disgnated pileup sites", false));

		// optionlist.add(getOption("cb", "chrBands", true,
		// "Chromosome Band",false));

		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	public void exec(CommandLine cl) throws Exception {
		// String mappability =
		// "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";

		List<String> files = new ArrayList<String>();
		String id = cl.getOptionValue("id");
		String normalbamf = cl.getOptionValue("n");
		String tumorbamf = cl.getOptionValue("t");
		String twobitref = cl.getOptionValue("r");
		String dbSNP = cl.getOptionValue("snp");
		String targetRegion = null;
		if (cl.hasOption("ct")) {
			targetRegion = cl.getOptionValue("ct");
		}
		String outdir = cl.getOptionValue("o");
		if (!outdir.endsWith("/")) {
			outdir = outdir + "/";
			File f = new File(outdir);
			if (!f.exists()) {
				boolean suc = false;
				try {
					suc = f.mkdirs();
				} catch (Exception ex) {
					System.out.println("could not make directory " + outdir);
					return;
				}
				if (suc == false) {
					System.out.println("could not make directory " + outdir);
					return;
				}
			}
		}
		//
		files.add(normalbamf);
		files.add(tumorbamf);
		files.add(twobitref);
		files.add(dbSNP);
		if(targetRegion!=null){
			files.add(targetRegion);
		}
		// files.add(outdir);

		boolean nopdf = cl.hasOption("nopdf");

		String mappability = null;
		if (cl.hasOption("mp")) {
			mappability = cl.getOptionValue("mp");
			files.add(mappability);
		}
		String prop = null;
		if (cl.hasOption("prop")) {
			prop = cl.getOptionValue("prop");
			files.add(prop);
		}
		String g1000 = null;
		if (cl.hasOption("g1000")) {
			g1000 = cl.getOptionValue("g1000");
			files.add(g1000);
		}
		String cosmic = null;
		if (cl.hasOption("cosmic")) {
			cosmic = cl.getOptionValue("cosmic");
			files.add(cosmic);
		}

		String exonSNP = null;
		if (cl.hasOption("exonSNP")) {
			exonSNP = cl.getOptionValue("exonSNP");
			files.add(exonSNP);
		}

		float g1000thres = 0f;
		if (cl.hasOption("g1000freq")) {
			g1000thres = Float.parseFloat(cl.getOptionValue("g1000freq"));
		}

		String readsStat = null;
		if (cl.hasOption("rs")) {
			readsStat = cl.getOptionValue("rs");
			String[] sa = readsStat.split(",");
			files.add(sa[0]);
			files.add(sa[1]);
		}

		String refflat = null;
		if (cl.hasOption("rg")) {
			refflat = cl.getOptionValue("rg");
			files.add(refflat);
		}

		String band = null;
		if (cl.hasOption("cb")) {
			band = cl.getOptionValue("cb");
			files.add(band);
		}

		String targetChr = null;
		boolean fullanalysis = true;
		if (cl.hasOption("chr")) {
			fullanalysis = false;
			targetChr = cl.getOptionValue("chr");
		}

		String startend = null;
		if (cl.hasOption("startend")) {
			startend = cl.getOptionValue("startend");
		}

		boolean useAvearageNormal = false;
		if (cl.hasOption("nd")) {
			useAvearageNormal = cl.getOptionValue("ng")
					.equalsIgnoreCase("true");
		}

		String sites = null;
		if (cl.hasOption("sites")) {
			sites = cl.getOptionValue("sites");
		}

		boolean allfileexsist = fileExsistanceCheck(files);
		if (!allfileexsist) {
			return;
		}
		// load property from file
		KarkinosProp.load(prop);
		String outobjdir = outdir + "sobj/";
		File fd = new File(outobjdir);
		fd.mkdirs();
		String outputsave = outobjdir + id + "saveobj.obj";

		SaveBean bean = null;
		File f = new File(outputsave);
		boolean loadFromSaved = f.exists();
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		String startends = startend;
		if(startends==null)startends="";

		if (targetChr != null) {
			if (!targetChr.contains("chr")) {
				String s = "chr" + targetChr;
				if (tgr.isRefExsist(s)) {
					targetChr = s;
				} else {
					String idxChr = targetChr;
					if (targetChr.contains("chr")) {
						idxChr = targetChr.replace("chr", "");
					}
					try {
						int chridxn = Integer.parseInt(idxChr);
						SAMFileReader normalbamr = getReader(normalbamf);
						targetChr = normalbamr.getFileHeader()
								.getSequence(chridxn - 1).getSequenceName();
						normalbamr.close();
					} catch (Exception ex) {
					}
					// try{
					// int chridxn = Integer.parseInt(idxChr);
					// targetChr = tgr.getChromString(chridxn);
					// }catch(Exception ex){
					// }
				}
			}
			outputsave = outobjdir + targetChr + "_" + id + "_"+ startends
					+ "_saveobj.obj";
		}

		if (loadFromSaved) {
			System.out.println("load data from pileuped middel file");
			bean = getReadsAndPileupDataFromFile(outputsave, tumorbamf, tgr);
		} else {
			System.out.println("load data from bam files");
			bean = getReadsAndPileupDataFromBam(normalbamf, tumorbamf, tgr,
					targetRegion, outputsave, targetChr,startend, refflat,sites);
		}
		// debug
		//fullanalysis=true;
		if (fullanalysis) {
			// analysis
			try {
				GeneExons ge = new GeneExons(refflat);
				// analysis(bean, dbSNP, mappability, tgr, tumorbamf, g1000,
				// g1000thres,cosmic, -1,exonSNP,ge);
				analysisNew(bean, dbSNP, mappability, tgr, normalbamf,tumorbamf, g1000,
						g1000thres, cosmic, -1, exonSNP, ge, useAvearageNormal,sites);
				// output
				output(bean, outdir, tgr, id, readsStat, alCNV, ge, na, pi,
						baseploidy, nopdf,refflat,sites);
			} catch (final IOException e) {
				System.out.println("could not read `" + refflat + "`: " + e.getMessage());
				return;
			}
		}
	}

	protected boolean fileExsistanceCheck(List<String> files) {
		for (String s : files) {
			File f = new File(s);
			if (!f.exists()) {
				System.out.println("file does not exsist " + s);
				return false;
			}
		}
		return true;
	}

	protected SaveBean getReadsAndPileupDataFromFile(String outputsave,
			String tumorbamf, TwoBitGenomeReader tgr) throws IOException,
			ClassNotFoundException {
		SAMFileReader tumorbamr = getReader(tumorbamf);
		List<SAMSequenceRecord> ssrList = tumorbamr.getFileHeader()
				.getSequenceDictionary().getSequences();
		tumorbamr.close();
		SaveBean sb = LoadSave.load(outputsave, ssrList, tgr);
		return sb;
	}

	AllelicCNV alCNV;
	NoiseAnalysis na;
	PeaksInfo pi;

	int baseploidy = 2;

	protected void analysisNew(SaveBean bean, String dbSNP, String mappability,
			TwoBitGenomeReader tgr,String normalbamf, String tumorbamf, String g1000,
			float g1000thres, String cosmic, double fixtumorratio,
			String exonSNP, GeneExons ge, boolean useAvearageNormal, String sites)
			throws IOException, ClassNotFoundException {
		System.out.println("analysis start");

		DataSet dataset = bean.getDataset();
		//debug

		dataset.setUseAvearageNormal(useAvearageNormal);
		ReadsSummary readsSummary = bean.getReadsSummary();

		// // Excute CNV analysis using wavelet transform
		// System.out.println("CNV analysis 1 start");
		GCParcentAdjust.calc(dataset);

		GeneEachCNV.calcCNVForEachgene(dataset,ge);

		denoizeWTAndHMM(dataset);
		// if CNV exceed max num 50 denoize more
		if (CountCNV.count(dataset) > 100) {
			KarkinosProp.maxdenoiseLevel = KarkinosProp.maxdenoiseLevel + 1;
			denoizeWTAndHMM(dataset);
		}
		// annotate datalist by CNV value
		// annotate datalist by CNV value
		dataset.assginCaptureInterval();

		int cnvcount = CountCNV.count(dataset);

		System.out.println("assgin dbSNP start");
		// // annotate datalist by dbSNP
		DbSNPAnnotation dbAnno = new DbSNPAnnotation(dbSNP, g1000, g1000thres,
				cosmic, exonSNP);
		bean.getDataset().assgindbSNP(dbAnno);
		// try to find allelic CNV
		alCNV = new AllelicCNV(dataset, readsSummary);
		//
		MatchMatrixBean matchmatrix = resolvePloidy(dataset, alCNV, pi,
				cnvcount);
		pi.setMatchmatrix(matchmatrix);
		double interval = matchmatrix.getInterval();
		if (interval > 0) {
			dataset.setBaselineLOHEstimate(interval);
		}

		FilterAnnotation fa = new FilterAnnotation(mappability, tgr, normalbamf,tumorbamf,
				dbAnno, ge);
		// //find hetro SNP correlation in each CNV interval
		// //reject if correration is high
		CorrelVaridate.varidate(dataset, alCNV);
		if (cnvcount < 300) {
			CheckPossibleHDAmp.check(dataset, pi, matchmatrix.getPloidyflg(),
					matchmatrix.getpEvenMax());
		}

		baseploidy = matchmatrix.getPloidyflg();
		CNVUtils.reflectToSNV(dataset, alCNV.getAllelicLOHLow(),
				alCNV.getAllelicLOHhigh());
		CorrelVaridate.recheck(dataset);

		dataset.getAnalyseDist();
		if (fixtumorratio > 0) {
			dataset.setFixtc((float) fixtumorratio);
		}

		// // support reads,entropy,mappability check
		System.out.println("filter annotation");

		// //set filter1
		float tc1 = dataset.getTumorRatio();
		fa.filterAnnotation(dataset, readsSummary, matchmatrix.getPloidyflg());
		// tcontents from somatic val
		dataset.getAnalyseDist().reanalyseTC(dataset);
		dataset.setBaseploidy(baseploidy);

		float tc2 = dataset.getTumorRatio();
		if (tc1 != tc2) {
			fa.filterAnnotation(dataset, readsSummary,
					matchmatrix.getPloidyflg());
		}
		// // set filter stat
		dataset.getAnalyseDist().analyseDist(dataset);

		dataset.bailAnalysis();
		System.out.println("analysis done");
		na = fa.getNa();
	}

	private MatchMatrixBean resolvePloidy(DataSet dataset, AllelicCNV alCNV2,
			PeaksInfo pi2, int cnvcount) {
		PloidyResolve pr = new PloidyResolve();
		pr.resolve(dataset, alCNV2, pi2, cnvcount);
		return pr.getMmb1();
	}

	private void denoizeWTAndHMM(DataSet dataset) throws IOException {
		System.out.println("CNV analysis 2 start");
		// MovingAverage.calc(dataset);
		double baselineLOHEstimate = WaveletDenoize.calc(dataset);
		// EMMethod to calc distribution
		System.out.println("CNV analysis 3 start");
		try{
			pi = EMMethod.calc(dataset, baselineLOHEstimate);
		}catch(Exception ex){
			ex.printStackTrace();
		}
		HMMCNVAnalysisFromEM.calc(dataset, pi);
	}

	protected void output(SaveBean bean, String outdir, TwoBitGenomeReader tgr,
			String id, String readsStat, AllelicCNV alCNV, GeneExons ge,
			NoiseAnalysis na2, PeaksInfo pi, int baseploidy, boolean nopdf,String refflat, String sites)
			throws Exception {
		DataSet dataset = bean.getDataset();
		float purity = dataset.getTumorRatio();
		ReadsSummary readsSummary = bean.getReadsSummary();
		// FileOutPut.outPutCNVData(outdir + id + "_cnvdata.txt", dataset);

		FileOutPut.lowcovBed(outdir + id + "_normal_lowcov.bed", outdir + id
				+ "_tumor_lowcov.bed", readsSummary,ge);

		TextSummary.outTextData(outdir + id + "_textdata.txt", readsSummary,
				readsStat, dataset, pi.getPloidy());

		CNVVcf.outData(outdir + id + "_cnvdata.vcf", dataset);

		FileOutPut.outPutSNVDataVCF(outdir + id + "_snvdata.vcf", dataset, tgr,
				na2);

		FileOutPut.outPutSNVDataForAnnover(outdir + id + "_annover_input.txt",
				dataset, tgr);

		FileOutPut.outputSNP(outdir + id + "_normalsnp.vcf", dataset, tgr, ge);

		FileOutPut.outputgeneCNV(outdir + id + "_cnv_forgene.txt", ge,refflat);

		if(sites!=null){
			FileOutPut.sites(outdir + id + "_disignatedsites.vcf", dataset, tgr, ge,sites);
		}

		FileOutPut.allDiff(outdir + id + "_alldiff.vcf", dataset, tgr, ge);

		if (!nopdf) {
			PdfReport.report(readsSummary, readsStat, dataset, alCNV, na2, pi,
					id, outdir + id + "_pdfdata.pdf");
		}
		try {
			FileOutPut.outputDepth(dataset, outdir + id + "_cnvdepth.txt",purity,ge);
		} catch (Exception ex) {
		}

		try {
			FileOutPut.outputAlleleDepth(alCNV, outdir + id
					+ "_cnvAllelicDepth.txt",purity);
		} catch (Exception ex) {
		}
	}

	public SaveBean getReadsAndPileupDataFromBam(String normalbamf,
			String tumorbamf, TwoBitGenomeReader tgr, String targetRegion,
			String outputsave, String targetChr, String startend,String refflat, String sites)
			throws Exception {
		DataSet dataset = new DataSet(tgr.readIndex());
		dataset.loadTargetBed(targetRegion, tgr);

		// Iterate by Chromosome
		SAMFileReader normalbamr = getReader(normalbamf);
		SAMFileReader tumorbamr = getReader(tumorbamf);

		ReadsSummary readsSummary = new ReadsSummary();
		readsSummary.setNormalbam(normalbamf);
		readsSummary.setTumorbam(tumorbamf);
		readsSummary.setRefFlat(refflat);
		readsSummary.setTaretbed(targetRegion);

		List<SAMSequenceRecord> ssrList = normalbamr.getFileHeader()
				.getSequenceDictionary().getSequences();
		int chromcnt = 0;

		System.out.println("start pileup " + targetChr);
		int cnt = 0;
		for (SAMSequenceRecord ssr : ssrList) {
			String chrom = ssr.getSequenceName();
			if (notEquals(targetChr, chrom)) {
				continue;
			}

			int length = ssr.getSequenceLength();
			if (tgr.isRefExsist(chrom)) {
				readsSummary.resetAlreadyregset();
				execChrom(dataset, chrom,startend, length, normalbamr, tumorbamr, tgr,
						readsSummary,sites);
				cnt++;
			} else {
				execChromOnlyForReadsStats(chrom, normalbamr, tumorbamr,
						readsSummary);
			}
			chromcnt++;
			readsSummary.resetAlreadyregset();
			System.gc();
		}
		dataset.getNormal().clearmap();
		dataset.getTumor().clearmap();

		if (cnt == 0) {
			System.out
					.println("no chromosome match bam to reference. check reference name");
		} else {
			System.out.println(cnt + " reference were processes");
		}

		// Basyan scoring
		// 1. collect LOH Hetero SNP

		// 2. collect 2N hetro

		// assgin by mappability

		// assgin neighbor indel and mutation check

		// output VCF, excel file

		// output graph

		// Save to Obj
		// no need to serialize
		readsSummary.clearGeneExon();
		SaveBean sbean = new SaveBean(dataset, readsSummary);
		LoadSave.save(sbean, outputsave);
		return sbean;
	}

	private boolean notEquals(String targetChr, String chrom) {
		if (targetChr == null) {
			return false;
		}
		String s1 = targetChr.replaceAll("chr", "");
		String s2 = chrom.replaceAll("chr", "");
		if (s1.equals(s2)) {
			return false;
		}
		return true;
	}

	private void execChromOnlyForReadsStats(String chrom,
			SAMFileReader normalbamr, SAMFileReader tumorbamr,
			ReadsSummary readsSummary) {
		CloseableIterator<SAMRecord> normarIte = null;

		try {
			normarIte = normalbamr.query(chrom, 0, 0, false);
		} catch (Exception ex) {
			return;
		}
		while (normarIte.hasNext()) {
			SAMRecord sam = normarIte.next();
			if (sam.getReadUnmappedFlag())
				continue;
			readsSummary.regN(sam, false, null);
		}
		normarIte.close();
		CloseableIterator<SAMRecord> tumorIte = tumorbamr.query(chrom, 0, 0,
				false);
		while (tumorIte.hasNext()) {
			SAMRecord sam = tumorIte.next();
			if (sam.getReadUnmappedFlag())
				continue;
			readsSummary.regT(sam, false, null);
		}
		tumorIte.close();
	}

	private void execChrom(DataSet dataset, String chrom, String startend, int length,
			SAMFileReader normalbamr, SAMFileReader tumorbamr,
			TwoBitGenomeReader tgr, ReadsSummary readsSummary, String sites)
			throws IOException {
		List<Interval> ivlist = ListUtils.getIntervalList(chrom, startend, length,
				KarkinosProp.BINBITSIZE);
		if(ivlist ==null){
			return;
		}
		DefinedSites ds = null;
		if(sites!=null){
			ds = new DefinedSites(sites);
			ds.load(sites,chrom);
		}
		//Ueda add 2020.12.18 load genome
		tgr.readRef(chrom);

		int n = 0;
		for (Interval iv : ivlist) {
			n++;
			//debug
			//if(n==10)break;

			CloseableIterator<SAMRecord> normarIte = normalbamr.query(chrom,
					iv.getStart(), iv.getEnd(), false);
			CloseableIterator<SAMRecord> tumorIte = tumorbamr.query(chrom,
					iv.getStart(), iv.getEnd(), false);

			int normalcnt = 0;
			int readslenn = 0;
			List<SamHolder> normalList = new ArrayList<SamHolder>();
			while (normarIte.hasNext()) {
				SAMRecord sam = normarIte.next();
				if (sam.getReadUnmappedFlag())
					continue;
				if (qualityCheck(sam)) {
					boolean onTarget = false;
					OntagetInfo oi = dataset.setNomalCoverageInfo(sam);
					if(oi!=null){
						onTarget = oi.ontag;
					}
					boolean dupli = sam.getDuplicateReadFlag();
					if (onTarget && !dupli) {
						SamHolder sh = new SamHolder();
						sh.setSam(sam);
						sh.setOi(oi);
						normalList.add(sh);
					}
					readsSummary.regN(sam, onTarget, iv);
				}
				normalcnt++;
				readslenn = sam.getReadLength();
			}

			System.out.println(iv.getStr() + " normal reads " + normalcnt
					+ " has been reads");

			List<SamHolder> tumorList = new ArrayList<SamHolder>();
			int tumorcnt = 0;
			int readslent = 0;
			while (tumorIte.hasNext()) {
				SAMRecord sam = tumorIte.next();
				//add 2020/12/17 for FFPE anneling near repeat, H.Ueda
				//extends softclip
				SoftClipExtention.extendSoftclip(sam, tgr);
				//count terminal mismatch
				if (sam.getIntegerAttribute("NM") != null && sam.getIntegerAttribute("NM") >= 2) {

					int terminalMismatch = TerminalMismatch.terminalMismatch(sam, tgr, KarkinosProp.extraReadTerminalCheckLen);
					if(terminalMismatch>=2){
						continue;
					}

				}
				if (qualityCheck(sam)) {
					boolean onTarget = false;
					OntagetInfo oi = dataset.setTumorCoverageInfo(sam);
					if(oi!=null){
						onTarget = oi.ontag;
					}
					boolean dupli = sam.getDuplicateReadFlag();
					if (onTarget && !dupli) {

						SamHolder sh = new SamHolder();
						sh.setSam(sam);
						sh.setOi(oi);
						tumorList.add(sh);
					}
					readsSummary.regT(sam, onTarget, iv);
				}
				tumorcnt++;
				readslent = sam.getReadLength();
			}

			System.out.println(iv.getStr() + " tumor reads " + tumorcnt
					+ " has been reads");

			dataset.setBinReads(iv, normalcnt, tumorcnt);
			System.out.println(Calendar.getInstance().getTime());
			PileUP.pileup(iv, dataset, normalList, tumorList, tgr, readsSummary,ds);
			System.out.println(Calendar.getInstance().getTime());

			System.out.println("normal depth"+ readsSummary.getNormalDepth().getMeanOntagDepth());
			System.out.println("tumor depth" + readsSummary.getTumorDepth().getMeanOntagDepth());

			readsSummary.setReadslent(readslent);
			readsSummary.setReadslenn(readslenn);
			normalList = null;
			tumorList = null;
			normarIte.close();
			tumorIte.close();
			System.gc();
		}
	}

	private boolean qualityCheck(SAMRecord sam) {
		return !SamUtils.lowmap(sam);
	}
}
