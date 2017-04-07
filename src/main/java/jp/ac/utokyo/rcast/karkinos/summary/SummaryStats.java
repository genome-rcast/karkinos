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
package jp.ac.utokyo.rcast.karkinos.summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import au.com.bytecode.opencsv.CSVReader;

public class SummaryStats {

	private static List<Option> getOption() {
		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("d", "directory", true,
				"input karkinos directory"
						+ " for each sample, cmmma separated", true));
		optionlist.add(getOption("o", "out", true, "output directory", true));
		optionlist.add(getOption("id", "id", true, "id for sampleset", false));
		optionlist.add(getOption("cb", "chromosomeBand", true,
				"chromosome band file from ucsc", false));

		optionlist.add(getOption("gr", "generef", true,
				"refflat base gene annotation files", false));

		optionlist.add(getOption("ct", "captureTarget", true,
				"Capture target regions(bed format)", false));

		optionlist
				.add(getOption(
						"rc",
						"minreccurent",
						true,
						"calculate p-val down to min reccurent sample per gene",
						false));

		optionlist.add(getOption("hg", "hg.2bit", true,
				"human genome .2bit reference", false));

		// optionlist.add(getOption("kg", "kgXref", true,
		// "kgXref file from ucsc", false));

		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	public static void main(String[] arg) throws IOException {

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

		String indir = cl.getOptionValue("d");
		String outdir = cl.getOptionValue("o");
		String id = "all";
		if (cl.hasOption("id")) {
			id = cl.getOptionValue("id");
		}

		String cb = null;
		if (cl.hasOption("cb")) {
			cb = cl.getOptionValue("cb");
		}

		String gr = null;
		if (cl.hasOption("gr")) {
			gr = cl.getOptionValue("gr");
		}

		String hg = null;
		if (cl.hasOption("hg")) {
			hg = cl.getOptionValue("hg");
		}

		String ct = null;
		long caplength = 60326149;
		if (cl.hasOption("ct")) {
			ct = cl.getOptionValue("ct");
			caplength = GetCparegionLength
					.loadTargetFromCapregionFile(new File(ct));
		}
		int minrecurrent = 2;
		if (cl.hasOption("rc")) {
			String minr = cl.getOptionValue("rc");
			minrecurrent = Integer.parseInt(minr);
		}

		List<Filebean> list = new ArrayList<Filebean>();

		for (String s : indir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive(f, list);
		}

		File dir = new File(outdir);
		if (!dir.exists()) {
			dir.mkdirs();
		}
		File outcvs = new File(dir.getAbsoluteFile() + "/" + id
				+ "_mutation_summary.csv");

		// write to cvs
		try {
			outToCSV(list, outdir, outcvs, id);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		File outexcel2 = new File(dir.getAbsoluteFile() + "/" + id
				+ "_after_final_filter_summary.xlsx");
		writeexcel(cb, list, outcvs, outexcel2, true, gr, caplength,
				minrecurrent, hg);

		File outexcel1 = new File(dir.getAbsoluteFile() + "/" + id
				+ "_summary.xlsx");
		writeexcel(cb, list, outcvs, outexcel1, false, gr, caplength,
				minrecurrent, hg);

		// analyse CNV by gene

		// write signature info
		// File substetutionInfo = new File(dir.getAbsoluteFile() + "/" + id
		// + "_substitutionBase.csv");
		//

	}

	private static boolean searchRecursive(File f, final List<Filebean> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("textdata.txt")
						|| name.contains("annoplus.csv");
				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive(new File(absPath), list);
				}
			}

		});

		if (listFiles == null)
			return false;

		List<File> flist = new ArrayList<File>();
		for (File ff : listFiles) {
			if (ff.isFile()) {
				flist.add(ff);
			}
		}

		for (File ff : flist) {

			if (ff.getName().contains("textdata.txt")) {

				Filebean fb = new Filebean();
				String sid = ff.getName().substring(0,
						ff.getName().indexOf("_textdata"));
				fb.setId(sid);
				fb.setTdata(ff);
				for (File ff2 : flist) {
					if (ff2.getName().contains("annoplus.csv")
							&& ff2.getName().startsWith(sid)) {

						fb.setCsvdata(ff2);

					}
				}
				if (fb.getCsvdata() != null) {
					System.out.println(fb.getTdata().getAbsolutePath());
					System.out.println(fb.getCsvdata().getAbsolutePath());
					System.out.println(fb.getId());
					list.add(fb);
				}
			}
		}
		return true;

	}

	private static void writeexcel(String cb, List<Filebean> list, File outcvs,
			File outexcel1, boolean bfilter, String gr, long caplength,
			int minrecurrent, String hg) {
		// register to memory DB
		try {

			XSSFWorkbook wb = new XSSFWorkbook();
			CellStyle cs1 = getCS(wb, 3);
			CellStyle cs2 = getCS(wb, 2);
			CellStyle cs3 = getCS(wb, 1);
			CellStyle[] csa = new CellStyle[] { cs1, cs2, cs3 };

			SummaryDB sb = new SummaryDB(outcvs);

			if (gr != null) {
				sb.loadRef(new File(gr), "generef");
			}

			DataReader dr = new DataReader(list);

			XSSFSheet statsheet = wb.createSheet("AA_cahnge_mutation_by_gene");
			sb.geneStat(dr, statsheet, bfilter, csa, true, wb, caplength,
					false, minrecurrent,false);

			XSSFSheet statsheetamphd = wb
					.createSheet("AA_cahnge_mutation_by_gene(AmpHD)");
			sb.geneStat(dr, statsheetamphd, bfilter, csa, true, wb, caplength,
					true, minrecurrent,false);

			XSSFSheet statsheet_all = wb.createSheet("all_mutation_by_gene");
			sb.geneStat(dr, statsheet_all, bfilter, csa, false, wb, caplength,
					true, minrecurrent,false);

			//
			XSSFSheet mutationstats = wb.createSheet("mutation_type");
			int colcnt = sb.mutationStat(mutationstats, 0, true, bfilter);
			sb.mutationStat(mutationstats, colcnt + 5, false, bfilter);

			if (hg != null) {
				XSSFSheet mutationsig = wb.createSheet("mutation_signature");				
				sb.mutationSig(mutationsig,hg);
			}
			//
			XSSFSheet tumorcontents = wb.createSheet("tumor_contents_ratio");

			DataWriter.writeTr(dr, tumorcontents, sb.getSampleL());

			ChromBand cband = null;
			if (cb != null) {
				cband = new ChromBand(cb);
			}

			XSSFSheet statsheet_HD_amp = wb.createSheet("HD_amp");
			sb.writeHDAMP(dr, statsheet_HD_amp, cband);

			//
			XSSFSheet cnvlist = wb.createSheet("cnvlist");
			DataWriter.writeCNV(dr, cnvlist, sb.getSampleL(), 1);

			XSSFSheet cnvlistal = wb.createSheet("cnvlistAllelic");
			DataWriter.writeCNV(dr, cnvlistal, sb.getSampleL(), 2);

			if (cb != null) {
				//
				XSSFSheet cnvcblist = wb.createSheet("cnvchrbandlist");
				DataWriter.writeCNVCB(dr, cnvcblist, sb.getSampleL(), cband,
						csa, false);

				XSSFSheet cnvcblistlow = wb.createSheet("cnvchrbandlist(low)");
				DataWriter.writeCNVCBAL(dr, cnvcblistlow, sb.getSampleL(),
						cband, csa, false);

				XSSFSheet cnvcblisthigh = wb
						.createSheet("cnvchrbandlist(high)");
				DataWriter.writeCNVCBAL(dr, cnvcblisthigh, sb.getSampleL(),
						cband, csa, true);

				XSSFSheet cnvcblistf = wb.createSheet("focal_cnvchrbandlist");
				DataWriter.writeCNVCB(dr, cnvcblistf, sb.getSampleL(), cband,
						csa, true);

				cband.close();
			}

			XSSFSheet rstat = wb.createSheet("readstats");
			DataWriter.writeReadsStats(dr, rstat, sb.getSampleL());

			sb.close();

			FileOutputStream out = new FileOutputStream(outexcel1);
			wb.write(out);
			out.close();

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static CellStyle getCS(XSSFWorkbook wb, int i) {
		CellStyle cs = wb.createCellStyle();
		cs.setFillPattern(CellStyle.SOLID_FOREGROUND);
		if (i == 1) {
			//
			cs.setFillForegroundColor(IndexedColors.LIGHT_YELLOW.getIndex());

		} else if (i == 2) {
			//
			cs.setFillForegroundColor(IndexedColors.GREEN.getIndex());

		} else {
			cs.setFillForegroundColor(IndexedColors.PINK.getIndex());
		}
		return cs;
	}

	private static void outToCSV(List<Filebean> list, String outdir,
			File outcvs, String id) throws IOException {

		File dir = new File(outdir);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outcvs));
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(
				dir.getAbsoluteFile() + "/" + id
						+ "_mutation_AAchange_summary.csv"));
		boolean isFirst = true;
		String[] csvtitle = null;
		int filteridx = 0;
		int FuncIdx = 0;
		int exonFuncIdx = 0;
		int geneIndex = 0;
		for (Filebean fb : list) {

			CSVReader brcvs = new CSVReader(new FileReader(fb.getCsvdata()));
			csvtitle = brcvs.readNext();

			if (isFirst) {
				isFirst = false;
				filteridx = index(csvtitle, "Filter");
				FuncIdx = index(csvtitle, "Func");
				exonFuncIdx = index(csvtitle, "ExonicFunc");
				geneIndex = index(csvtitle, "Gene");
				removeSpace(csvtitle);
				writeLine(bw, csvtitle, -1, csvtitle.length);
				writeLine(bw2, csvtitle, -1, csvtitle.length);

			}

			String[] data = null;
			while ((data = brcvs.readNext()) != null) {

				if (filterOK(data, filteridx)) {

					writeLine(bw, data, geneIndex, csvtitle.length);
					if (aaChange(data, FuncIdx, exonFuncIdx)) {
						writeLine(bw2, data, geneIndex, csvtitle.length);
					}

				}

			}

		}
		bw.close();
		bw2.close();

	}

	private static void removeSpace(String[] csvtitle) {
		int n = 0;
		for (String s : csvtitle) {
			s = s.replaceAll(" ", "");
			csvtitle[n] = s;
			n++;
		}

	}

	private static int index(String[] data, String s2) {

		int n = 0;
		for (String s : data) {

			if (s.equalsIgnoreCase(s2)) {
				return n;
			}
			n++;
		}
		return n;
	}

	private static boolean filterOK(String[] data, int filteridx) {

		String filter = data[filteridx];
		boolean filOK = filter.equalsIgnoreCase("PASS");
		return filOK;

	}

	public static boolean aaChange(String[] data, int func, int exonFuncIdx) {
		if (data.length - 1 < exonFuncIdx) {
			return true;
		}
		String funcS = data[func].trim();
		String efunc = data[exonFuncIdx].trim();
		//
		boolean AAChange = funcS.contains("exonic")
				|| funcS.contains("splicing");
		if (funcS.contains("ncRNA")) {
			AAChange = false;
		}
		if ((efunc.equalsIgnoreCase("synonymous SNV"))) {
			AAChange = false;
		}
		// boolean noAAChange = (efunc.length() < 2) ||
		// (efunc.equalsIgnoreCase("synonymous SNV"));
		// boolean AAChange = !noAAChange;
		return AAChange;
	}

	private static void writeLine(BufferedWriter bw, String[] ary,
			int geneIndex, int length) throws IOException {

		StringBuffer sb = new StringBuffer();
		int n = 0;
		if (ary != null) {
			for (String s : ary) {
				s = s.replaceAll(",", ";");
				if (n == geneIndex) {
					s = "\"" + s + "\"";
				}
				if (sb.length() == 0) {
					sb.append(s);
				} else {
					sb.append("," + s);
				}
				n++;
				if (n == length - 1)
					break;
			}
		}
		bw.write(sb.toString() + "\n");

	}

}
