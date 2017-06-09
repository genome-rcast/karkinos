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
package jp.ac.utokyo.rcast.karkinos.exec.develop;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.karkinos.ploidy.PloidyResolve;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.LoadSave;
import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.SaveBean;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.graph.output.CNVVcf;
import jp.ac.utokyo.rcast.karkinos.graph.output.FileOutPut;
import jp.ac.utokyo.rcast.karkinos.graph.output.PdfReport;
import jp.ac.utokyo.rcast.karkinos.graph.output.TextSummary;
import jp.ac.utokyo.rcast.karkinos.hmm.CountCNV;
import jp.ac.utokyo.rcast.karkinos.hmm.HMMCNVAnalysis;
import jp.ac.utokyo.rcast.karkinos.hmm.HMMCNVAnalysisFromEM;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
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

public class TestTG extends ReadWriteBase {

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
			help.printHelp("karkinos.jar reanalysis", opts, true);
			return;
		}

		TestTG tg = new TestTG();
		try {
			tg.exec(cl);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	private static List<Option> getOptionListForKarkinos() {

		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("md", "middledate", true,
				"middle data object", true));
		optionlist
				.add(getOption("t", "tumorBam", true, "tumor bam file", true));
		optionlist.add(getOption("r", "reference", true,
				"2 bit genome reference file", true));
		optionlist.add(getOption("snp", "dbSNP", true,
				"dbSNP list from annover sites,(bin,chr,start,end)", true));
		optionlist.add(getOption("ct", "captureTarget", true,
				"Capture target regions(bed format)", true));

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

		optionlist.add(getOption("rs", "readsStats", true,
				"reads stats files(normal,tumor)", false));

		optionlist.add(getOption("rg", "refFlatGenes", true,
				"optional,gene reference for depth stats", false));

		optionlist.add(getOption("tc", "tumorContents", true,
				"fixed tumor contents (and skip tc calculation)", false));

		optionlist.add(getOption("exonSNP", "exonSNP", true,
				"additional Exon SNP", false));

		// optionlist.add(getOption("cb", "chrBands", true,
		// "Chromosome Band",false));

		return optionlist;
	}

	public void exec(CommandLine cl) throws Exception {

		// String mappability =
		// "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";

		System.out.println("Analysis starts with");
		System.out.println(KarkinosProp.getInfoString());
		List<String> files = new ArrayList<String>();
		String id = cl.getOptionValue("id");
		String middledata = cl.getOptionValue("md");
		String tumorbamf = cl.getOptionValue("t");
		String twobitref = cl.getOptionValue("r");
		String dbSNP = cl.getOptionValue("snp");
		String targetRegion = cl.getOptionValue("ct");
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
		files.add(middledata);
		files.add(tumorbamf);
		files.add(twobitref);
		files.add(dbSNP);
		files.add(targetRegion);
		// files.add(outdir);

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
		String band = null;
		if (cl.hasOption("cb")) {
			band = cl.getOptionValue("cb");
			files.add(band);
		}

		String refflat = null;
		if (cl.hasOption("rg")) {
			refflat = cl.getOptionValue("rg");
			files.add(refflat);
		}

		String tc = null;
		double dtc = -1;
		if (cl.hasOption("tc")) {
			tc = cl.getOptionValue("tc");
			try {
				dtc = Double.parseDouble(tc);
				if (dtc > 1 || dtc < 0.1) {
					System.out
							.println("tumor contents ratio have to be between 0.1 to 1");
					return;
				}
			} catch (Exception ex) {
				System.out
						.println("tumor contents ratio have to be between 0.1 to 1");
				return;
			}
		}

		boolean allfileexsist = fileExsistanceCheck(files);
		if (!allfileexsist) {
			return;
		}
		// load property from file
		KarkinosProp.load(prop);

		SaveBean bean = null;
		File f = new File(middledata);
		boolean loadFromSaved = f.exists();
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		 if (loadFromSaved) {
		 System.out.println("load date from pileuped middel file");
		 bean = getReadsAndPileupDataFromFile(middledata, tumorbamf);
		
		 //
		
		
		 }
		GeneExons ge = new GeneExons(refflat);
		// analysis
		analysisNew(bean, dbSNP, mappability, tgr, tumorbamf, g1000,
				g1000thres, cosmic, dtc, exonSNP, ge);
		// output
		// output(bean, outdir, tgr, id, readsStat, alCNV,ge,na,pi);

	}

	AllelicCNV alCNV;
	NoiseAnalysis na;
	PeaksInfo pi;

	protected void analysisNew(SaveBean bean, String dbSNP, String mappability,
			TwoBitGenomeReader tgr, String tumorbamf, String g1000,
			float g1000thres, String cosmic, double fixtumorratio, String exonSNP,GeneExons ge) throws IOException, ClassNotFoundException {

		System.out.println("analysis start");
		DataSet dataset = bean.getDataset();
		ReadsSummary readsSummary = bean.getReadsSummary();
		

		GCParcentAdjust.calc(dataset);
		denoizeWTAndHMM(dataset);
		//if CNV exceed max num 25 denoize more 
		if(CountCNV.count(dataset)>25){
			KarkinosProp.maxdenoiseLevel = KarkinosProp.maxdenoiseLevel+1;
			denoizeWTAndHMM(dataset);
		}
		// annotate datalist by CNV value
		dataset.assginCaptureInterval();		
		System.out.println("assgin dbSNP start");
//		// annotate datalist by dbSNP
//		DbSNPAnnotation dbAnno = new DbSNPAnnotation(dbSNP, g1000, g1000thres,cosmic,exonSNP);
//		bean.getDataset().assgindbSNP(dbAnno);
	    //save save

//		save(bean,data1);
//		save(dbAnno,data2);



	}

	private void resolvePloidy(DataSet dataset, AllelicCNV alCNV2, PeaksInfo pi2) {
		//
		PloidyResolve pr = new PloidyResolve();
		pr.resolve(dataset, alCNV2, pi2,20);

	}

	public static Object _load(File f) throws IOException,
			ClassNotFoundException {

		// Read from disk using FileInputStream
		FileInputStream f_in = new FileInputStream(f);

		// Read object using ObjectInputStream
		ObjectInputStream obj_in = new ObjectInputStream(f_in);

		// Read an object
		Object obj = obj_in.readObject();
		return obj;

	}

	public static void save(Object sbean, String outputsave) throws IOException {

		// Write to disk with FileOutputStream
		FileOutputStream f_out = new FileOutputStream(outputsave);

		// Write object with ObjectOutputStream
		ObjectOutputStream obj_out = new ObjectOutputStream(f_out);

		// Write object out to disk
		obj_out.writeObject(sbean);

	}

	protected void output(SaveBean bean, String outdir, TwoBitGenomeReader tgr,
			String id, String readsStat, AllelicCNV alCNV, GeneExons ge,
			NoiseAnalysis na2, PeaksInfo pi) throws Exception {

		DataSet dataset = bean.getDataset();
		ReadsSummary readsSummary = bean.getReadsSummary();
		// FileOutPut.outPutCNVData(outdir + id + "_cnvdata.txt", dataset);

		FileOutPut.lowcovBed(outdir + id + "_normal_lowcov.bed", outdir + id
				+ "_tumor_lowcov.bed", readsSummary,null);

		TextSummary.outTextData(outdir + id + "_textdata.txt", readsSummary,
				readsStat, dataset,pi.getPloidy());
		CNVVcf.outData(outdir + id + "_cnvdata.vcf", dataset);
		FileOutPut.outPutSNVDataVCF(outdir + id + "_snvdata.vcf", dataset, tgr,
				na2);
		FileOutPut.outPutSNVDataForAnnover(outdir + id + "_annover_input.txt",
				dataset, tgr);

		FileOutPut.outputSNP(outdir + id + "_normalsnp.vcf", dataset, tgr, ge);

		PdfReport.report(readsSummary, readsStat, dataset, alCNV, na2, pi,id,
				outdir + id + "_pdfdata.pdf");

	}

	private void _denoizeAndHMM(DataSet dataset) throws IOException {
		System.out.println("CNV analysis 2 start");
		WaveletDenoize.calc(dataset);

		System.out.println("CNV analysis 3 start");
		// HMM call for CNV
		HMMCNVAnalysis.calc(dataset);
	}

	// for future use
	private void denoizeAndHMM(DataSet dataset) throws IOException {
		System.out.println("CNV analysis 2 start");
		// MovingAverage.calc(dataset);
		double baselineLOHEstimate = WaveletDenoize.calc(dataset);
		// EMMethod to calc distribution
		System.out.println("CNV analysis 3 start");
		pi = EMMethod.calc(dataset, baselineLOHEstimate);
		System.out.println("CNV analysis 4 start");
		// HMM call for CNV
		HMMCNVAnalysis.calc(dataset);

	}

	private void denoizeWTAndHMM(DataSet dataset) throws IOException {

		System.out.println("CNV analysis 2 start");
		// MovingAverage.calc(dataset);
		double baselineLOHEstimate = WaveletDenoize.calc(dataset);
		// EMMethod to calc distribution
		System.out.println("CNV analysis 3 start");
		pi = EMMethod.calc(dataset, baselineLOHEstimate);
		HMMCNVAnalysisFromEM.calc(dataset, pi);

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
			String tumorbamf) throws IOException, ClassNotFoundException {

		SAMFileReader tumorbamr = getReader(tumorbamf);
		List<SAMSequenceRecord> ssrList = tumorbamr.getFileHeader()
				.getSequenceDictionary().getSequences();
		tumorbamr.close();
		SaveBean sb = LoadSave.load(outputsave, ssrList,null);
		return sb;

	}

}
