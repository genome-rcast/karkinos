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

public class MergeSample {

	public static void main(String[] arg) {

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

		Map<String,String[]> cvslines = new HashMap<String,String[]>();

		// reads cvs
		String[] csvtitle = null;

		try {
			
			//BufferedReader brcvs = new BufferedReader(new FileReader(incvs));
			CSVReader brcvs = new CSVReader(new FileReader(incvs));
			csvtitle = brcvs.readNext();
						
			String[] data =null; 
			//csvtitle = line.split(",");
			
			int index = 0;
			while ((data =brcvs.readNext())!= null) {
				

				
				String chr = data[getIndex(csvtitle,"Chr")].trim();
				//System.out.println(chr);
				String pos = data[getIndex(csvtitle,"Start")].trim();
				String end = data[getIndex(csvtitle,"End")].trim();
				String ref = data[getIndex(csvtitle,"Ref")].trim();
				String alt = data[getIndex(csvtitle,"Obs")].trim();


				int npos = 0;
				try{
					npos = Integer.parseInt(pos);
				}catch(NumberFormatException nfe){
					//
					//System.out.println(pos);
					
				}
				///
				boolean del = alt.equals("-");
				boolean ins = ref.equals("-");
				if(del){
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
			vcftitle = new String[] { 
					"chrom", "pos",
					"id", "REF", "ALT",
					"QUAL","Filter"};
			writeLine(bw, vcftitle, infotitle, csvtitle);
			while (line != null) {
				line = brvcf.readLine();
				if(line==null)break;
				if (line.startsWith("#")) {
					continue;
				}
				String[] data = line.split("\t");
				int pos = Integer.parseInt(data[1]);
				int len = data[4].length()-1;
				String key = data[0] + "-" + (pos+len);

				String[] vfcdata = data;
				String infoline = vfcdata[7];
				String[] vfcdataw = new String[vfcdata.length-1];
				int idx = 0;
				for(String s:vfcdata){
					if(idx==vfcdata.length-1)break;
					vfcdataw[idx]=vfcdata[idx];
					idx++;
				}
				String[] infodata = getInfoStr(infoline);
				String[] cvsdata = cvslines.get(key);
				if(cvsdata==null){
					key = data[0] + "-" + (pos);
					cvsdata = cvslines.get(key);
				}
				if(cvsdata==null){
					key = data[0] + "-" + (pos-1);
					cvsdata = cvslines.get(key);
				}
				
				writeLine(bw, vfcdataw, infodata, cvsdata);
			}
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static int getIndex(String[] title, String string) {
		
		int idx = 0;
		for(String s:title){
			
			if(s.equalsIgnoreCase(string)){
				return idx;
			}
			idx++;
			
		}
		return idx;
	}

	static final String[] infotitle = new String[] { 
			"Depth",
			"Allele Freq",
			"Allel Freq Org", 
			"CopyNumber",
			"RMS BaseQuality",
			"RMS MapQuality",
			"Fisher p val", 
			"Seq Entropy",
			"Mappability",
			"Fisher P val for direction Check",
			"LOGn",
			"LOGt",
			"LOGt Adjuated",
			"dbSNP",
			"Repeat",
			"Supported by One direction" };

	private static String[] getInfoStr(String infoline) {
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
			Map<String,String> map = new HashMap<String,String>();
			for(String s:dev2){
				String[] sa = s.split("=");
				map.put(sa[0],sa[1]);				
			}
			
			ret[0] = map.get("DP");
			ret[1] = map.get("AF");
			ret[2] = map.get("AFO");
			ret[3] = map.get("CN");
			ret[4] = map.get("BQ");
			ret[5] = map.get("MQ");
			ret[6] = map.get("pval");
			ret[7] = map.get("S");
			ret[8] = map.get("MP");
			ret[9] = map.get("pD");
			ret[10] = map.get("logn");
			ret[11] = map.get("logt");
			ret[12] = map.get("logtAdjuated");

			ret[13] = infoset1.contains("DB") ? "Yes" : "";
			ret[14] = infoset1.contains("RP") ? "Yes" : "";
			ret[15] = infoset1.contains("OD") ? "Yes" : "";
			
		} catch (Exception ex) {

		}
		return ret;
	}

	private static void writeLine(BufferedWriter bw, String[] ary0,
			String[] ary1, String[] ary2) throws IOException {

		StringBuffer sb = new StringBuffer();
		if (ary0 != null) {
			for (String s : ary0) {
				s=s.replaceAll(",", ";");
				sb.append(s + ",");
			}
		}
		if (ary1 != null) {
			for (String s : ary1) {
				if(s!=null){
					s=s.replaceAll(",", ";");
				}
				sb.append(s + ",");
			}
		}
		if (ary2 != null) {
			for (String s : ary2) {
				if(s!=null){
					s=s.replaceAll(",", ";");
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
		return optionlist;
	}

	public static Option getOption(String opt, String longOpt, boolean hasArg,
			String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

}
