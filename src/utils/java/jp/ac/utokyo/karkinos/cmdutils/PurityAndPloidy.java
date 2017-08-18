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
import java.util.ArrayList;
import java.util.List;

public class PurityAndPloidy {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_NCC/karkinosver4.0.25_2,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor/,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/tcga,"
				+ "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/goss,";

		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/pulityploidy.txt";
		FileWriter fw = new FileWriter(new File(out));
		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter pw = new PrintWriter(bw);
		
		
		List<File> list = new ArrayList<File>();
		for (String s : indir.split(",")) {
			File f = new File(s);
			searchRecursive(f, list);
		}
		for (File f : list) {

			//
			exec(f,pw);
			
		}
		pw.close();
		
		
	}

	private static void exec(File f, PrintWriter pw) throws IOException {
		
		String id = f.getName();
		id = id.replaceAll("_textdata.txt", "");
		//
		float purity = 0f;
		float ploidy = 0f;
		
		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(f));
		while ((s = br.readLine()) != null) {
			
			String[] sa = s.split("\t");
			if(s.startsWith("tc used")){
				purity = Float.parseFloat(sa[1]);
			}
			if(s.startsWith("ploidy")){
				ploidy = Float.parseFloat(sa[1]);
			}
			
		}
		System.out.println(id+"\t"+purity+"\t"+ploidy);
		pw.print(id+"\t"+purity+"\t"+ploidy);
		pw.print("\n");
		
		
	}

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("_textdata.txt");
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
			if (ff.isFile() && ff.getName().contains("_textdata.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

}
