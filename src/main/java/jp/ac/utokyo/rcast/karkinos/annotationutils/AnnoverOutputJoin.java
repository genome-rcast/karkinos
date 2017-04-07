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
package jp.ac.utokyo.rcast.karkinos.annotationutils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import au.com.bytecode.opencsv.CSVReader;

public class AnnoverOutputJoin {

	public static void main(String[] arg) {

		exec1(arg);
		// exec2(arg);

	}

	public static void exec1(String[] arg) {

		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOption();
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
			help.printHelp("karkinos.jar mergeresult", opts, true);
			return;
		}

		String incvs = cl.getOptionValue("a");
		String invcf = cl.getOptionValue("v");

		//
		File incvsf = new File(incvs);

		//
		String name = incvsf.getName();
		String name2 = name.substring(0, name.lastIndexOf('.'));
		String outpath = incvsf.getParentFile().getAbsolutePath() + "/" + name2
				+ "annoplus.csv";

		String id = name2;
		int testidx = id.indexOf("_annover");
		if (testidx > 0) {
			id = id.substring(0, testidx);
		}

		Map<String, String[]> cvslines = new HashMap<String, String[]>();

		// reads cvs
		String[] csvtitle = null;

		try {

			// BufferedReader brcvs = new BufferedReader(new FileReader(incvs));
			CSVReader brcvs = new CSVReader(new FileReader(incvs));
			csvtitle = brcvs.readNext();

			String[] data = null;
			// csvtitle = line.split(",");

			int index = 0;
			while ((data = brcvs.readNext()) != null) {

				String chr = data[getIndex(csvtitle, "Chr")].trim();
				// System.out.println(chr);
				String pos = data[getIndex(csvtitle, "Start")].trim();
				String end = data[getIndex(csvtitle, "End")].trim();
				String ref = data[getIndex(csvtitle, "Ref")].trim();
				String alt = data[getIndex(csvtitle, "Obs")].trim();

				int npos = 0;
				try {
					npos = Integer.parseInt(pos);
				} catch (NumberFormatException nfe) {
					//
					// System.out.println(pos);

				}
				// /
				boolean del = alt.equals("-");
				boolean ins = ref.equals("-");
				if (del) {
					npos = npos - 1;
				}

				String key = chr + "-" + npos;
				cvslines.put(key, data);

				index++;
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// reads vcf
		String[] vcftitle = null;
		try {

			BufferedWriter bw = new BufferedWriter(new FileWriter(outpath));
			BufferedReader brvcf = new BufferedReader(new FileReader(invcf));

			String line = brvcf.readLine();
			vcftitle = new String[] { "chrom", "pos", "id", "REF", "ALT",
					"QUAL", "Filter", "Filter2", "Normal Ref:Alt",
					"Tumor Ref:Alt", "pre seq", "post seq" };
			writeLine(bw, "sample_ID", vcftitle, infotitle, csvtitle);

			while (line != null) {
				line = brvcf.readLine();
				if (line == null)
					break;
				if (line.startsWith("#")) {
					continue;
				}
				String[] data = line.split("\t");
				int pos = Integer.parseInt(data[1]);
				int len = data[4].length() - 1;
				String key = data[0] + "-" + (pos + len);

				String[] vfcdata = data;

				String[] vfcdataw = new String[vcftitle.length];
				int idx = 0;
				for (String s : vfcdata) {
					if (idx == vcftitle.length - 1)
						break;
					vfcdataw[idx] = vfcdata[idx];
					idx++;
				}
				// for filter2
				vfcdataw[7] = vfcdata[8];
				if (vfcdata.length > 9) {
					vfcdataw[8] = vfcdata[9];
					vfcdataw[9] = vfcdata[10];
					if (vfcdata.length > 11) {
						vfcdataw[10] = vfcdata[11];
						vfcdataw[11] = vfcdata[12];
					}
				} else {
					vfcdataw[8] = "";
					vfcdataw[9] = "";
					vfcdataw[10] = "";
					vfcdataw[11] = "";
				}

				String infoline = vfcdata[7];
				String[] cvsdata = cvslines.get(key);
				int posn = pos;
				if (cvsdata == null) {
					key = data[0] + "-" + (pos);
					cvsdata = cvslines.get(key);
				}
				if (cvsdata == null) {
					key = data[0] + "-" + (pos - 1);
					cvsdata = cvslines.get(key);
					posn = pos - 1;
				}
				if (cvsdata == null) {
					key = data[0] + "-" + (pos + 1);
					cvsdata = cvslines.get(key);
					posn = pos + 1;
				}

				String geneSymbol = getGs(cvsdata);
				String[] infodata = getInfoStr(infoline, vfcdata[6],
						geneSymbol, data[0], posn);
				writeLine(bw, id, vfcdataw, infodata, cvsdata);

			}
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static String getGs(String[] cvsdata) {
		if (cvsdata == null || cvsdata.length < 2) {
			return "";
		} else {

			String gs = cvsdata[1];
			if (gs.length() < 2) {
				System.out.println(gs);
			}
			if (gs.indexOf('(') > 0) {
				gs = gs.substring(0, gs.indexOf('('));
			}
			if (gs.indexOf(';') > 0) {
				gs = gs.substring(0, gs.indexOf(';'));
			}
			return gs;

		}

	}

	public static void exec2(String[] arg) {

		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOption();
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
			help.printHelp("karkinos.jar mergeresult", opts, true);
			return;
		}

		String incvs = cl.getOptionValue("a");
		String invcf = cl.getOptionValue("v");
		File incvsf = new File(incvs);
		String name = incvsf.getName();
		String name2 = name.substring(0, name.lastIndexOf('.'));

		String inpath = incvsf.getParentFile().getAbsolutePath() + "/" + name2
				+ "annoplus.csv";

		String outpath = incvsf.getParentFile().getAbsolutePath() + "/" + name2
				+ "_final.maf";

		// /

	}

	private static int getIndex(String[] title, String string) {

		int idx = 0;
		for (String s : title) {

			if (s.equalsIgnoreCase(string)) {
				return idx;
			}
			idx++;

		}
		return idx;
	}

	static final String[] infotitle = new String[] { "SNV pval", "chr_a",
			"pos_a", "geneSymbol", "highprob_SNP", "initfail", "cosmic",
			"cosmic_validated", "cosmic_count", "Indel", "Depth",
			"Allele Freq", "Allele Freq Org", "CopyNumber",
			"A allele CopyNumber", "B Allele CopyNumber", "RMS BaseQuality",
			"RMS MapQuality", "Fisher p val", "Seq Entropy", "Mappability",
			"Fisher P val for direction Check", "LOGn", "LOGt",
			"LOGn Adjuated", "LOGt Adjuated", "Normal Depth", "Normal Ratio",
			"Second Mutation Allele", "Second mutation AF", "dbSNP",
			"dbSNP validation status", "Repeat", "Supported by One direction",
			"Allelic info Available", "reference B allele frequlecy",
			"mutation B allle frequency","score", };

	private static String[] getInfoStr(String infoline, String passFileters,
			String geneSymbol, String chr, int pos) {
		String[] ret = new String[infotitle.length];
		try {
			Set<String> infoset1 = new HashSet<String>();
			String[] dev1 = infoline.split(";");
			if (dev1.length > 1) {
				int n = 0;
				for (String s : dev1) {
					if (n == 0) {
						n++;
						continue;
					}
					infoset1.add(s);
					n++;
				}
			}
			String infos2 = dev1[0];
			String[] dev2 = infos2.split(",");
			Map<String, String> map = new HashMap<String, String>();
			for (String s : dev2) {
				if (s.contains("=")) {
					String[] sa = s.split("=");
					if (sa.length == 2) {
						map.put(sa[0], sa[1]);
					}

				}
			}

			// String filter2 = passFileters;
			// if(infoset1.contains("LowRefFilter")){
			// filter2 = "lowref";
			// }
			// if(infoset1.contains("LowTumorFilter")){
			// if(filter2.equals("PASS")){
			// filter2 = "lowtumor";
			// }else if(filter2.length()==0){
			// filter2 = "lowtumor";
			// }else{
			// filter2 = filter2 +",lowtumor";
			// }
			// }
			String bq = map.get("BQ");
			boolean indel = false;
			if (bq != null) {
				bq = bq.trim();
				float bqf = 0f;
				try {
					bqf = Float.parseFloat(bq);
					indel = bqf == 0;
				} catch (Exception ex) {
				}
			}
			boolean cosmic = infoset1.contains("cosmic");
			boolean cosmic_v = infoset1.contains("cosm_v");
			boolean initfail = infoset1.contains("initfail");
			boolean validateSNP = infoset1.contains("validateSNP");
			ret[0] = map.get("pvSNV");
			ret[1] = chr;
			ret[2] = String.valueOf(pos);
			ret[3] = geneSymbol;
			ret[4] = validateSNP ? "Yes" : "";
			ret[5] = initfail ? "Yes" : "";
			ret[6] = cosmic ? "Yes" : "";
			ret[7] = cosmic_v ? "Yes" : "";
			ret[8] = nullempty(map.get("cosmc"));

			ret[9] = indel ? "Yes" : "";
			ret[10] = map.get("DP");
			ret[11] = map.get("AF");
			ret[12] = map.get("AFO");
			ret[13] = map.get("CN");
			ret[14] = map.get("ACN");
			ret[15] = map.get("BCN");
			ret[16] = map.get("BQ");
			ret[17] = map.get("MQ");
			ret[18] = map.get("pval");
			ret[19] = map.get("S");
			ret[20] = map.get("MP");
			ret[21] = map.get("pD");
			ret[22] = map.get("logn");
			ret[23] = map.get("logt");
			ret[24] = map.get("lognAdjuated");
			ret[25] = map.get("logtAdjuated");
			ret[26] = map.get("ND");
			ret[27] = map.get("NR");

			ret[28] = map.get("SA");
			ret[29] = map.get("SAF");

			ret[30] = infoset1.contains("DB") ? "Yes" : "";
			ret[31] = nullempty(map.get("DBv"));
			ret[32] = infoset1.contains("RP") ? "Yes" : "";
			ret[33] = infoset1.contains("OD") ? "Yes" : "";

			ret[34] = infoset1.contains("AIA") ? "Yes" : "";
			ret[35] = nullempty(map.get("sBAF"));
			ret[36] = nullempty(map.get("mBAF"));

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return ret;
	}

	private static String nullempty(String s) {
		if (s == null) {
			return "";
		}
		return s;
	}

	private static void writeLine(BufferedWriter bw, String string,
			String[] ary0, String[] ary1, String[] ary2) throws IOException {

		StringBuffer sb = new StringBuffer();
		sb.append(string + ",");
		if (ary0 != null) {
			for (String s : ary0) {
				s = s.replaceAll(",", ";");
				sb.append(s + ",");
			}
		}
		if (ary1 != null) {
			for (String s : ary1) {
				if (s != null) {
					s = s.replaceAll(",", ";");
				}
				sb.append(s + ",");
			}
		}
		if (ary2 != null) {
			for (String s : ary2) {
				if (s != null) {
					s = s.replaceAll(",", ";");
				}
				sb.append(s + ",");
			}
		}
		//
		bw.write(sb.toString() + "\n");

	}

	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("v", "vcf", true, "vcf output of karkinos",
				true));
		optionlist.add(getOption("a", "annoverSummary", true,
				"annover summary output", true));
		optionlist.add(getOption("id", "sample_ID", true, "sample ID", false));
		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

}
