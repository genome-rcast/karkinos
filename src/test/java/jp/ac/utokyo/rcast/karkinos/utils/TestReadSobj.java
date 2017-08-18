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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map.Entry;

import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.LoadUtils;
import jp.ac.utokyo.rcast.karkinos.annotation.loadsave.SaveBean;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;

public class TestReadSobj {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String middelfile = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/challenge3/sobj";
		try {
			load(middelfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static SaveBean load(String s) throws IOException, ClassNotFoundException {

		File f = new File(s);

		//
		if (f.isDirectory()) {
			//
			List<File> flist = new ArrayList<File>();
			for (File fa : f.listFiles()) {
				if (fa.getName().endsWith(".obj")) {
					flist.add(fa);
				}
			}

			SaveBean sb = new SaveBean(null, null);
			sort(flist);
			int cnt = 0;
			int normaltotal = 0;
			int tumortotal = 0;
			
			for (File ff : flist) {

				 if(! (ff.getName().contains("chr1_")||ff.getName().contains("chr2_")) ){
					 //debug
					 continue;
				 }

				cnt++;
				System.out.println("reading " + ff.getName());

				SaveBean sbeach = null;
				try {
					sbeach = _load(ff);
				} catch (Exception ex) {

				}
				if (sbeach != null) {

					System.out.println("normal="
							+ sbeach.getReadsSummary().getNormalCounter()
									.getTotalmap());
					System.out.println("tumor="
							+ sbeach.getReadsSummary().getTumorCounter()
									.getTotalmap());
					
					normaltotal = normaltotal + sbeach.getReadsSummary().getNormalCounter()
					.getTotalmap();
					tumortotal = tumortotal + sbeach.getReadsSummary().getTumorCounter()
					.getTotalmap();
					
				}
				if(sbeach!=null){
					LoadUtils.merge(sb, sbeach);
				}

				// debug
				// for(SNVHolder snv:sbeach.dataset.getSnvlist()){
				// if(snv.getPos()==1418474){
				// System.out.println(snv.getChr()+"\t"+snv.getPos()+"\t"+snv.getFlg());
				// }
				// }
				// List<SNVHolder> list = sbeach.dataset.getSnvlist();
				// if (list != null && list.size() > 0) {
				// String chr = list.get(0).getChr();
				// if
				// (tgr.isRefExsist(chr)||tgr.isRefExsist(chr.replaceAll("chr","")))
				// {
				// System.out.println(chr);

	
				//

			}
			
			for (Entry<String, ReadsCounter> entry : sb.getReadsSummary().getNormalPerChrom().entrySet()) {

				String chrom = entry.getKey();
				ReadsCounter rc = entry.getValue();
				int target = rc.getTotalOnTarget();
				int total = rc.getTotalmap();
				int outoftarget = total - target;
				int totalu = rc.getTotalunique();
				int uniqueoutoftarget = totalu - rc.getTotalUniqueOntarget();

				System.out.println(chrom);
				System.out.println(target);
				System.out.println(total);
				

			}
			
			
			

		}
		return null;

	}
	
	private static void sort(List<File> flist){
		
		Collections.sort(flist,new Mycomp2());

	}

	static class Mycomp2 implements Comparator<File>{

		public int compare(File f0, File f1) {
			
			String s0 = f0.getName().substring(0, 4);
			String s1 = f1.getName().substring(0, 4);
			if(s0.equals(s1)){
				int n = getN(f0);
				int m = getN(f1);	
				return n-m;
				
			}else{
				return s0.compareTo(s1);
			}
			
			
		}

		private int getN(File f) {
			
			String s =  f.getName().substring(f.getName().indexOf("challenge3_")+10,f.getName().lastIndexOf("_"));
			s = s.replaceAll("_","");
			return Integer.parseInt(s);
		}
		
	}
	
	public static SaveBean _load(File f) throws IOException,
			ClassNotFoundException {

		// Read from disk using FileInputStream
		FileInputStream f_in = new FileInputStream(f);

		// Read object using ObjectInputStream
		ObjectInputStream obj_in = new ObjectInputStream(f_in);

		// Read an object
		Object obj = obj_in.readObject();
		SaveBean saveBean = null;
		if (obj instanceof SaveBean) {
			saveBean = (SaveBean) obj;
		}
		return saveBean;

	}
	
}
