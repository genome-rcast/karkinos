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
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.summary.CNAInterval;
import jp.ac.utokyo.rcast.karkinos.summary.ChromBand;

public class CompWithArray {

	public static void main(String[] arg) throws IOException, SQLException{
		
		//
		String cbhg19f = "/GLUSTER_DIST/data/users/ueda/genome/hg19chrband.txt";				
		String cbhg18f = "/GLUSTER_DIST/data/users/ueda/genome/hg18chrband.txt";	
		ChromBand cbhg18 = new ChromBand(cbhg18f,"hg18");
		ChromBand cbhg19 = new ChromBand(cbhg19f,"hg19");
		String out = "/GLUSTER_DIST/data/users/ueda/tesCNVtout.txt";
		//
		
		String array  = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/SNP6/selected/RCAST/txt/";
		String exome  = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast";
		
		List<File> list = new ArrayList<File>();
		searchRecursive(new File(exome), list);
		
		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);
		
		for(File f:list){
			
			int n = f.getName().indexOf("-");
			String id = f.getName().substring(0,n);
			File arrayf = getFile(array,id);
			if(arrayf==null){
				System.out.println(id+"\t"+f.getName()+"\t"+"null");
				continue;
			}else{
			//
				System.out.println(id+"\t"+f.getName()+"\t"+arrayf.getName());
				//
				exec(id,f,arrayf,cbhg18,cbhg19,pw);
			
			}
		}
		//
		pw.close();
		
		
	}
	
	private static void exec(String id, File f, File arrayf, ChromBand cbhg18,
			ChromBand cbhg19, PrintWriter pw) throws IOException, SQLException {
		
		Map<String, Bean> mapArray = new LinkedHashMap<String, Bean>();
		Map<String, Bean> mapExom = new LinkedHashMap<String, Bean>();
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String s = null;
		while ((s = br.readLine()) != null) {
			
			if(s.startsWith("#"))continue;
			String[] sa = s.split("\t");
			float highs = Float.parseFloat(sa[5]);
			float lows = Float.parseFloat(sa[6]);
			if (highs > 4)
				highs = 4;
			if (lows > 4)
				lows = 4;
			Bean b = null;
			String chr0 = sa[0];
			if (!chr0.contains("chr")) {
				chr0 = "chr" + chr0;
			}
			int pos = Integer.parseInt(sa[1]);
			String key = chr0 + cbhg19.getBand(chr0, pos, pos);
			if (mapExom.containsKey(key)) {
				b = mapExom.get(key);
			} else {
				b = new Bean();
				mapExom.put(key, b);
			}
			b.high.addValue(highs);
			b.low.addValue(lows);
			
		}

		BufferedReader abr = new BufferedReader(new FileReader(arrayf));
		String as = null;
		while ((as = abr.readLine()) != null) {
			String[] asa = as.split("\t");
			float highs = Float.parseFloat(asa[5]);
			float lows = Float.parseFloat(asa[6]);
			if (highs > 4)
				highs = 4;
			if (lows > 4)
				lows = 4;
			Bean b = null;
			String chr0 = "chr"+asa[1];
			int pos = Integer.parseInt(asa[2]);
			String key = chr0 + cbhg18.getBand(chr0, pos, pos);
			if (mapArray.containsKey(key)) {
				b = mapArray.get(key);
			} else {
				b = new Bean();
				mapArray.put(key, b);
			}
			b.high.addValue(highs);
			b.low.addValue(lows);
		}
		
		
		for(CNAInterval cna:cbhg19.getList()){
			
			String id0 = cna.getChr()+cna.getName();
			Bean arrayB = mapArray.get(id0);	
			Bean exomeB = mapExom.get(id0);	
			float highA = 1;
			float lowA =1;
			float highE =1;
			float lowE =1;
			if(arrayB!=null&&arrayB.high.getN()>100){
				highA = (float)arrayB.high.getMean();
				lowA = (float)arrayB.low.getMean();
			}
			if(exomeB!=null&&exomeB.high.getN()>100){
				highE = (float)exomeB.high.getMean();
				lowE = (float)exomeB.low.getMean();
			}
			if((highA!=1)&&(lowA!=1)&&(highE!=1)&&(lowE!=1)){
							
				pw.print(id+"\t"+id0+"\t"+highA+"\t"+lowA+"\t"+highE+"\t"+lowE+"\n");
			
			}			
			
			
						
			
			
		}
		
		
		
	}

	private static File getFile(String array,String id) {
		
		File f = new File(array);
		for(File ff:f.listFiles()){
			if(ff.getName().startsWith(id)){
				return ff;
			}
		}		
		return null;
	}

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("cnvAllelicDepth.txt");
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
		for (File ff : listFiles) {
			if (ff.isFile() && ff.getName().contains("cnvAllelicDepth.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}
	
	
}
