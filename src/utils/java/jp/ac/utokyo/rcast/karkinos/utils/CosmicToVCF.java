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
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class CosmicToVCF {

	public static void main(String[] arg) {

		TwoBitGenomeReader tgr = new TwoBitGenomeReader(
				new File(
						"/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/hg38.2bit"));
		String targetRegion = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/S3035822_Covered_sm_lift_hg38_padded_sm.bed";

		DataSet dataset;

		String dir = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/CosmicCompleteTargetedScreensMutantExport.tsv";
		
		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/CosmicCompleteTargetedScreensMutantExport.vcf";

		
		
				
				
				
		
		try {

			
			File filew = new File(out);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filew)));

			
			pw.write("##fileformat=VCFv4.1 \n");
			pw.write("##reference=/files/grch38/cosmic/v80/CosmicCompleteTargetedScreensMutantExport.tsv.gz");
			pw.write("#CHROM POS     ID        REF    ALT     QUAL FILTER INFO");
			
			dataset = new DataSet(tgr.readIndex());
			dataset.loadTargetBed(targetRegion, tgr);

			InputStream is = null;
			Reader r = null;
			BufferedReader br = null;
			
			int cnt=0,writecnt=0;
			try {
				//is = new GZIPInputStream(new BufferedInputStream(new FileInputStream(dir)));
				is = new BufferedInputStream(new FileInputStream(dir));
				r  = new InputStreamReader(is, "MS932");
				br = new BufferedReader(r);
				int lc = 0;
				for (;;) {
					String text = br.readLine();	//改行コードは含まれない
					if (text == null) break;
					if(lc==0){
						lc++;
						continue;
					}
					
					System.out.println(text);
					
					String[] ary = text.split("\t");
					String info = ary[11];
					String genomeversion = ary[22];
					if(!genomeversion.equals("38")){
						continue;
					}					
					
					String pos = ary[23];
					
					String chrom = pos.split(":")[0];
					String poss = pos.split(":")[1].split("-")[0];
							
					String strand = ary[24]; 
					String refalt = ary[17]; 
					
					int n =  Integer.parseInt(poss);
					String ref = tgr.getGenomicSeq("chr"+chrom, n, n, true);
					String alt = refalt.substring(refalt.indexOf(">")+1);
					if(strand.equals("-")){
						if(alt.equals("A")){
							alt ="T";
						}else if(alt.equals("T")){
							alt ="A";
						}else if(alt.equals("G")){
							alt ="C";
						}else if(alt.equals("C")){
							alt ="G";
						}
					}
					String id = ary[1]+"."+lc;
					String line = (chrom+"\t"+ poss+"\t"+id+"\t"+ref+"\t"+alt+"\t"+"100 \t PASS \t "+info);
					System.out.println(line);
					pw.write(line+"\n");
					lc++;
					
					//if(lc==100)break;
					
					
//					if(text.startsWith("#")){
//						
//					
//						pw.write(text+"\n");
//						
//						continue;
//						
//					}
//					String sa[] = text.split("\t");
//					String chr = "chr"+sa[0];
//					int pos = Integer.parseInt(sa[1]);
//					//
//					cnt++;
//					if(dataset.getCh().getCapInterval(chr, pos)!=null){
//						
//						pw.write(text+"\n");
//						 writecnt++;
//												
//					}
//					
//
//					if(cnt%1000==0){
//						System.out.println(writecnt +"/" +cnt +" " + chr +" "+pos);
//						
//					}						
					
				}
			} catch (Exception e) {
				//throw new RuntimeException(e);
				e.printStackTrace();
			} finally {
				if (br != null) try { br.close(); } catch (IOException e) {}
				if (r  != null) try { r .close(); } catch (IOException e) {}
				if (is != null) try { is.close(); } catch (IOException e) {}
			}
		

			
				

			pw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
