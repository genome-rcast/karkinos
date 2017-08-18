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

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class RefFlatFilter {

	public static void main(String[] arg) {

		String in = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/refFlat.txt.gz";

		String list = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/nm_ids.txt";

		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/refFlat_panel.txt";

		Set<String> s = new HashSet<String>();

		try {

			File file = new File(list);
			BufferedReader br0 = new BufferedReader(new FileReader(file));

			String str = br0.readLine();
			while (str != null) {
				System.out.println(str);
				s.add(str);
				str = br0.readLine();
			}

			br0.close();

			File filew = new File(out);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(
					filew)));

			InputStream is = null;
			Reader r = null;
			BufferedReader br = null;

			int cnt = 0, writecnt = 0;

			Map<String,List<String[]>> holder = new HashMap<String,List<String[]>>();
			
			String title ="";
			try {
				is = new GZIPInputStream(new BufferedInputStream(
						new FileInputStream(in)));
				// is = new BufferedInputStream(new FileInputStream(in));
				r = new InputStreamReader(is, "MS932");
				br = new BufferedReader(r);
				int lc = 0;
				for (;;) {
					String text = br.readLine(); // 改行コードは含まれない
					if (text == null)
						break;
					if (text.startsWith("#")) {
						
						title = text;
						continue;
					}

					String ar[] = text.split("\t");

					if (s.contains(ar[1])) {
						//pw.write(text + "\n");
						//reg
						String gene = ar[0];
						if(holder.containsKey(gene)){
							holder.get(gene).add(ar);
						}else{
							List<String[]> sl = new ArrayList<String[]>();
							sl.add(ar);
							holder.put(gene, sl);							
						}
						
					}
					lc++;

					// if(lc==100)break;

					// if(text.startsWith("#")){
					//
					//
					// pw.write(text+"\n");
					//
					// continue;
					//
					// }
					// String sa[] = text.split("\t");
					// String chr = "chr"+sa[0];
					// int pos = Integer.parseInt(sa[1]);
					// //
					// cnt++;
					// if(dataset.getCh().getCapInterval(chr, pos)!=null){
					//
					// pw.write(text+"\n");
					// writecnt++;
					//
					// }
					//
					//
					// if(cnt%1000==0){
					// System.out.println(writecnt +"/" +cnt +" " + chr
					// +" "+pos);
					//
					// }

				}
				
				if(title.length()>10){
					pw.write(title+"\n");
				}
				//
				Iterator<String> ite = holder.keySet().iterator();
				//
				while(ite.hasNext()){
					
					//
					String key = ite.next();
					List<String[]> ll = holder.get(key);
					String[] gcd = gcd(ll); 
					
					int n= 0;
					for(String sss:gcd){
						
						if(n>0){
							pw.write("\t");
						}
						pw.write(sss);
						n++;
						
					}
					pw.write("\n");
					
				}
				
				
			} catch (Exception e) {
				// throw new RuntimeException(e);
				e.printStackTrace();
			} finally {
				if (br != null)
					try {
						br.close();
					} catch (IOException e) {
					}
				if (r != null)
					try {
						r.close();
					} catch (IOException e) {
					}
				if (is != null)
					try {
						is.close();
					} catch (IOException e) {
					}
			}

			pw.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private static String[] gcd(List<String[]> ll) {
		
		if(ll.size()==1)return ll.get(0);
		
		String[] ret = new String[ll.get(0).length]; 
		//
		StringBuffer ids = new StringBuffer();
		int maxexon = 0;
		String exons = "";
		String exone = "";
		
		int txs = 0;
		int txe =0;
		int cdss = 0;
		int cdse = 0;
		
		int cnt=0;
		for(String[] ar:ll){
			
			if(cnt==0){
				
				txs = Integer.parseInt(ar[4]);
				txe = Integer.parseInt(ar[5]);
				cdss = Integer.parseInt(ar[6]);
				cdse = Integer.parseInt(ar[7]);
				
			}
			
			ret = ar;
			ids.append(ar[1]+",");
			
			int maxe = Integer.parseInt(ar[8]);
			if(maxe>maxexon){
				maxexon = maxe;
				exons = ar[9];
				exone = ar[10];
			}
			
			if(Integer.parseInt(ar[4])<txs)txs =Integer.parseInt(ar[4]);
			if(Integer.parseInt(ar[5])>txe)txe =Integer.parseInt(ar[5]);
			if(Integer.parseInt(ar[6])<cdss)cdss =Integer.parseInt(ar[6]);
			if(Integer.parseInt(ar[7])>cdse)cdse =Integer.parseInt(ar[7]);
			
			
			cnt++;
			
		}
		
		ret[1] = ids.toString();
		ret[4] = txs+"";
		ret[5] = txe+"";
		ret[6] = cdss+"";
		ret[7] = cdse+"";
		
		ret[8] = maxexon+"";
		ret[9] = exons;
		ret[10] = exone;
		
		//
		return ret;
	}

}
