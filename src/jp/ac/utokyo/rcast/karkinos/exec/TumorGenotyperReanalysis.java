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

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.SaveBean;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.utils.DataHolder;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class TumorGenotyperReanalysis extends TumorGenotyper {

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

		TumorGenotyperReanalysis tg = new TumorGenotyperReanalysis();
		try {
			tg.exec(cl);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private static List<Option> getOptionListForKarkinos() {

		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("md", "middledate", true,
				"middle data object", true));

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

		optionlist.add(getOption("rs", "readsStats", true,
				"reads stats files(normal,tumor)", false));

		optionlist.add(getOption("rg", "refFlatGenes", true,
				"optional,gene reference for depth stats", false));

		optionlist.add(getOption("tc", "tumorContents", true,
				"fixed tumor contents (and skip tc calculation)", false));

		optionlist.add(getOption("exonSNP", "exonSNP", true,
				"additional Exon SNP", false));

		optionlist.add(getOption("nd", "normaldepth", true,
				"use averaged normal depth if available ", false));

		// optionlist.add(getOption("cb", "chrBands", true,
		// "Chromosome Band",false));

		optionlist.add(getOption("nopdf", "nopdf", false,
				"no graphic summary pdf output", false));
		

		optionlist.add(getOption("sites", "pileupsites", true,
				"disgnated pileup sites", false));


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
		files.add(middledata);
		files.add(tumorbamf);
		files.add(twobitref);
		files.add(dbSNP);
		if (targetRegion != null) {
			files.add(targetRegion);
		}
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
		boolean nopdf = cl.hasOption("nopdf");

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
		boolean useAvearageNormal = false;
		if (cl.hasOption("nd")) {
			useAvearageNormal = cl.getOptionValue("nd")
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

		SaveBean bean = null;
		File f = new File(middledata);
		boolean loadFromSaved = f.exists();
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		if (loadFromSaved) {
			System.out.println("load date from pileuped middel file");
			bean = getReadsAndPileupDataFromFile(middledata, tumorbamf, tgr);
		}
		GeneExons ge = new GeneExons(refflat);
		// analysis
		// analysis(bean, dbSNP, mappability, tgr, tumorbamf, g1000,
		// g1000thres,cosmic,
		// dtc,exonSNP,ge);
		analysisNew(bean, dbSNP, mappability, tgr, normalbamf,tumorbamf, g1000,
				g1000thres, cosmic, dtc, exonSNP, ge, useAvearageNormal,sites);

		// output
		output(bean, outdir, tgr, id, readsStat, alCNV, ge, na, pi, baseploidy,
				nopdf,refflat,sites);

	}

}
