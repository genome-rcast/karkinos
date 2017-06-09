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
package jp.ac.utokyo.rcast.karkinos.annotation.loadsave;

import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

//to separate bam pileup and reads counting process to
//later analysis process;
public class LoadSave {

	public static void save(SaveBean sbean, String outputsave)
			throws IOException {

		// Write to disk with FileOutputStream
		FileOutputStream f_out = new FileOutputStream(outputsave);

		// Write object with ObjectOutputStream
		ObjectOutputStream obj_out = new ObjectOutputStream(f_out);

		// Write object out to disk
		obj_out.writeObject(sbean);

	}

	public static SaveBean load(String s, List<SAMSequenceRecord> ssrList,
			TwoBitGenomeReader tgr) throws IOException, ClassNotFoundException {

		File f = new File(s);
		if (f.isFile()) {
			return _load(s);
		}
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
			sort(flist, ssrList);
			int cnt = 0;
			for (File ff : flist) {

//				if(!ff.getName().contains("chr1_")){
//					//debug
//					continue;
//				}
				
				cnt++;				
				System.out.println("reading " + ff.getName());
				
				SaveBean sbeach = null;
				try{
					sbeach = _load(ff);
				}catch(Exception ex){
					
				}
				if(sbeach != null){
					
					System.out.println("normal="+sbeach.getReadsSummary().getNormalCounter().getTotalmap());
					System.out.println("tumor="+sbeach.getReadsSummary().getTumorCounter().getTotalmap());
				}
				
				// debug
				// for(SNVHolder snv:sbeach.dataset.getSnvlist()){
				// if(snv.getPos()==1418474){
				// System.out.println(snv.getChr()+"\t"+snv.getPos()+"\t"+snv.getFlg());
				// }
				// }
				//List<SNVHolder> list = sbeach.dataset.getSnvlist();
				//if (list != null && list.size() > 0) {
					//String chr = list.get(0).getChr();
					//if (tgr.isRefExsist(chr)||tgr.isRefExsist(chr.replaceAll("chr",""))) {
					//	System.out.println(chr);
					
				if(sbeach!=null){
					LoadUtils.merge(sb, sbeach);
				}
					//}
				//}
				//

			}
			return sb;

		}
		return null;

	}

	private static void sort(List<File> flist, List<SAMSequenceRecord> ssrList) {

		Collections.sort(flist, new FileComparator(ssrList));

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

	public static SaveBean _load(String f) throws IOException,
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
