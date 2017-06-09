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

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;



public class ReadsStatsVirus extends ReadWriteBase {
	
	
	public static void main(String[] arg){
		
		List<File> list = new ArrayList<File>();
		String bamfiledir = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_Baylor_1";
		for (String s : bamfiledir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive(f, list);
		}
		//
		
		
		for(File f: list){
			
			System.out.println(f.getName());
			SAMFileReader bamr = getReader(f);
			for(SAMSequenceRecord ssr:bamr.getFileHeader().getSequenceDictionary().getSequences()){
			
				
				if(ssr.getSequenceName().contains("chr")){
					continue;
				}else{
					int cnt=0;
					CloseableIterator<SAMRecord> ite = bamr.query(ssr.getSequenceName(),
							0, 0, false);
					
					while(ite.hasNext()){
						ite.next();
						cnt++;
					}					
					ite.close();
					if(ssr.getSequenceName().contains("HBV")){
						System.out.println(ssr.getSequenceName()+"\t"+cnt);
					}
				}				
				
			}
			
			
			

			
			
		}
		
		
	}
	
	
	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("tumor_genome.bam");
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
			if (ff.isFile() && ff.getName().contains("tumor_genome.bam")) {

				// System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

	
}
