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
package jp.ac.utokyo.karkinos.cmdutils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;
import java.util.TreeSet;

import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;

public class GermSNPToAnnover {

	public static void main(String[] arg) throws IOException {

		//
		String in = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNP.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNPAnnov.txt";

		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);

		//
		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(in));
		Set<String> set = new TreeSet<String>();

		while ((s = br.readLine()) != null) {

			String[] sa = s.split("\t");
			if(sa[7].contains("rs"))continue;
			String ref = sa[8];
			if(ref.contains(",")){
				ref = ref.substring(0,ref.indexOf(","));
			}
			
			String alt = sa[9];
			if(alt.contains(",")){
				alt = alt.substring(0,alt.indexOf(","));
			}
			
			String posid = sa[5] + "-" + sa[6] + "-" + ref + "-" + alt;
			set.add(posid);

		}
		br.close();

		for (String ss : set) {

			// System.out.println(ss);
			String[] sa = ss.split("-");
			//
			boolean del = sa[2].length() > 1;
			boolean ins = sa[3].length() > 1;

			if (del) {
							
				
				int start = Integer.parseInt(sa[1]);
				int end = start + (sa[2].length()-1);
				String line = sa[0] + "\t" + start + "\t" + end + "\t" + sa[2] + "\t"
				+ sa[3];
				pw.println(line);
				
				
			} else {

				String line = sa[0] + "\t" + sa[1] + "\t" + sa[1] + "\t" + sa[2] + "\t"
						+ sa[3];
				pw.println(line);
			}
		}

		pw.close();

	}

}
