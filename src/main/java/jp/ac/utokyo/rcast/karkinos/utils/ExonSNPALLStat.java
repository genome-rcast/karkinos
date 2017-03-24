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
package jp.ac.utokyo.rcast.karkinos.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class ExonSNPALLStat {

	public static void main(String[] arg) throws IOException{
		
		//
		Set<String> s = new HashSet<String>();
		
		String fs0 = "/GLUSTER_DIST/data/users/ueda/project/sakata/datar.txt";
		FileInputStream fis0 = new FileInputStream(new File(fs0));
		BufferedReader br0 = new BufferedReader(new InputStreamReader(fis0));
		for (;;) {
			
			String line = br0.readLine();
			if (line == null)
				break;
			if (line.startsWith("#")) {
				continue;
			}
			//
			s.add(line);
			
		}	
		
		String fs = "/GLUSTER_DIST/data/users/ueda/ICGC/exomeSNP5.vcf";
		FileInputStream fis = new FileInputStream(new File(fs));
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));
		//
		Map<Integer,int[]> cm = new TreeMap<Integer,int[]>();
		int cn= 0;
		for (;;) {
			
			String line = br.readLine();
			if (line == null)
				break;
			if (line.startsWith("#")) {
				continue;
			}
			String[] sa = line.split("\t");
			int n = Integer.parseInt(sa[5]);
			int[] ia = null;
			if(cm.containsKey(n)){
				ia = cm.get(n);
			}else{
				ia = new int[2];
				cm.put(n, ia);
			}
			if(sa[2].contains("rs")){
				ia[0] = ia[0]+1;	
			}else{
				ia[1]=ia[1]+1;
			}
			//
			cn++;
			int pos = Integer.parseInt(sa[1]);
			boolean disp = false;
			String pid = "chr"+sa[0]+"_"+sa[1];
			if(s.contains(pid)){
				System.out.println(line);
			}			
			
			
		}	
		
//		Iterator<Integer> ite = cm.keySet().iterator();
//		//
//		while(ite.hasNext()){
//			
//			//
//			int n= ite.next();
//			System.out.println(n+"\t"+cm.get(n)[0]+"\t"+cm.get(n)[1]);
//			
//		}
		
	}
}
