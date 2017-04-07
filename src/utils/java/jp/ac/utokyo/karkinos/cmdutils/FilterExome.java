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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class FilterExome {
	
	
	public static void main(String[] arg) throws IOException{
		
		String in1 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNP.txt";
		String in2 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNPAnnov.txt.exome_summary.csv";
		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/germlineSNPaanopuls.txt";
		
		//
		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);
		//
		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(in2));
		Map<String,String> map = new HashMap<String,String>();
		while ((s = br.readLine()) != null) {
					
			String[] sa = s.split(",");
			//System.out.println(sa[21]+"\t"+sa[22]+"\t"+sa[23]);	
			String posid = sa[21] + "-" + sa[22] + "-" + sa[24] + "-" + sa[25];			
			//System.out.println(posid);
			map.put(posid,s);
		
		}	
		
		br = new BufferedReader(new FileReader(in1));
		while ((s = br.readLine()) != null) {
					
			String[] sa = s.split("\t");
			String ref = sa[8];
			if(ref.contains(",")){
				ref = ref.substring(0,ref.indexOf(","));
			}
			if(sa[7].contains("rs"))continue;
			
			String alt = sa[9];
			if(alt.contains(",")){
				alt = alt.substring(0,alt.indexOf(","));
			}
			
			String posid = sa[5] + "-" + sa[6] + "-" + ref + "-" + alt;
			//
			if(map.containsKey(posid)){
				//
				pw.println(s+"\t"+map.get(posid));
				
			}			
		
		}	
		
		
		pw.close();
		
	}
	

}
