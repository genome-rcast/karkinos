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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class TestSummaryStat {

	/**
	 * @param args 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
	
		
		//String s1 = "/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_baylor,/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_rcast";
		
		//String s1 = "/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_baylor,/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_rcast,/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_nccall";
		
//		String s1 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast," +
//				"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_NCC/karkinosver4.0.25," +
//				"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor";
		
		String s1 = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT";
		
//		String s1 = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/TN," 
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT16-1-MT16N,"
//		+  	"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT16-2-MT16N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT22-1-MT22N,"
//		+  	"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT22-2-MT22N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT36-1-MT36N,"
//		+  	"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT36-2-MT36N,"
//		+  	"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT10_2-TN10N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT34-1-MT34N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT34-2-MT34N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT34-4-MT34N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT34-5-MT34N,"
//		+		"/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT34-6-MT34N";
		
//		String s1 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast," +
//		"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_NCC/karkinosver4," +
//		"/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor";
		
		//String s1 = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor";
		
		String s2 = "/GLUSTER_DIST/data/users/ueda/project/Aihara/20140530Exome/";
		String s3 = "MTExome";
		String s4 ="/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";
		
		String s5 ="/usr/local/karkinos/karkinos/genome/refgenes_anno.csv";
		
		String s6 ="/usr/local/karkinos/karkinos/genome/captureregionv4all.bed.capregion";
		
		String s7 ="/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		
		String recurrent = "2";
		
//		String s1 = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST,/GLUSTER_DIST/data/users/ueda/project/ICGC_HCC";
//		
//		String s2 = "/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/TESTSummary2";
//		String s3 = "TEST";
//		String s4 ="/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";
//		
//		String s5 ="/usr/local/karkinos/karkinos/genome/refgenes_anno.csv";
//		
//		String s6 ="/usr/local/karkinos/karkinos/genome/captureregionv4all.bed.capregion";
		
		
//		String s = "/GLUSTER_DIST/data/users/yamamoto/exome/Brain/MT";
//		File f = new File(s);
//		StringBuffer sb = new StringBuffer();
//		for(File ff:f.listFiles()){
//			
//			if(!ff.isDirectory())continue;
//			String ss = ff.getAbsolutePath();
//			if(sb.length()>0){
//				sb.append(",");
//			}
//			sb.append(ff.getAbsolutePath());
//			
//		}		
//		
//		String s1 = sb.toString();		
//		System.out.println(s1);
//		String s2 = "/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/MT_summary";
//		String s3 = "MT";
		
		//arg[3] = "-cb ";
		
		List<String> l = new ArrayList<String>();		
		add(l,"-d",s1);
		add(l,"-o",s2);
		add(l,"-id",s3);
		add(l,"-cb",s4);
		
		add(l,"-gr",s5);
		add(l,"-ct",s6);
		add(l,"-rc",recurrent);
		
		add(l,"-hg",s7);
		
		
		String[] ar = l.toArray(new String[l.size()]);
		SummaryStats.main(ar);
		//SummaryStatsVaridate.main(ar);

	}
	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

}
