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
package jp.ac.utokyo.rcast.karkinos.cmd;

import java.io.IOException;

import jp.ac.utokyo.rcast.karkinos.annotationutils.AnnoverOutputJoin;
import jp.ac.utokyo.rcast.karkinos.exec.TumorGenotyper;
import jp.ac.utokyo.rcast.karkinos.exec.TumorGenotyperReanalysis;
import jp.ac.utokyo.rcast.karkinos.summary.SummaryStats;
import jp.ac.utokyo.rcast.karkinos.summary.SummaryStatsVaridate;
import jp.ac.utokyo.rcast.karkinos.varidation.VaridateByBam;
//import jp.ac.utokyo.rcast.realign.Realignment;

public class KarkinosCmd {

	//public static final String version = "Karkinos2.0 version 0.1 2013/4/22";
	//public static final String version = "Karkinos2.0 version 0.3 2016/11/30";
	public static final String version = "Karkinos2.0 version 0.4 2017/03/06";

	public static void main(String[] arg) throws Exception {

		//
		// arg = new String[]{"analysis"};
		if (arg == null || arg.length == 0) {
			printMsg();
		} else {

			String cmd = arg[0];
			String[] arg2 = getArg(arg);
			if (cmd.equals("analysis")) {

				TumorGenotyper.main(arg2);

			} else if (cmd.equals("reanalysis")) {

				TumorGenotyperReanalysis.main(arg2);

			} else

			if (cmd.equals("mergeresult")) {

				AnnoverOutputJoin.main(arg2);

			} else if (cmd.equals("Summary")) {

				SummaryStats.main(arg2);
			
//			} else if (cmd.equals("realign")) {
//
//				Realignment.main(arg2);	
				
			} else if (cmd.equals("varidate")) {

				VaridateByBam.main(arg2);

			} else if (cmd.equals("SummaryVaridated")) {

				SummaryStatsVaridate.main(arg2); 
			
			}else {  

				printMsg();
			}

		}

	}

	private static String[] getArg(String[] arg) {
		String[] arg2 = new String[arg.length - 1];
		for (int n = 1; n < arg.length; n++) {
			arg2[n - 1] = arg[n];
			// System.out.println(arg[n]);
		}
		return arg2;
	}

	private static void printMsg() {
		System.out.println(version);
		System.out.println("usage: karkinos.jar <command> options");
		System.out.println("Command: ");
		System.out
				.println("		 analysis	analysis SNV,CSV,tumor rate from bamfiles");
		System.out
				.println("		 reanalysis	reanalysis of SNV,CSV,tumorrate from middle pileuped file");
		System.out
		.println("		 realign	realignment normar,tumor bam file");
		
		System.out.println("		 mergeresult	merge Annover result into csv file");
		System.out.println("		 Summary	summary result across sample");
		System.out.println("		 varidate	varidate SNV result by bam file");

	}

}
