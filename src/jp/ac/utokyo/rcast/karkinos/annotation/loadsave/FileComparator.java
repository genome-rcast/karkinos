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
import java.util.Comparator;
import java.util.List;

public class FileComparator implements Comparator<File> {

	List<SAMSequenceRecord> ssrList = null;

	public FileComparator(List<SAMSequenceRecord> _ssrList) {
		ssrList = _ssrList;
	}

	public int compare(File f1, File f2) {
		int idx1 = getIndex(ssrList, f1);
		int idx2 = getIndex(ssrList, f2);
		if (idx1 != idx2) {
			return idx1 - idx2;
		} else {
			int n = comp(f1,f2);
			return n;
		}
	}

	private int getIndex(List<SAMSequenceRecord> ssrList2, File f1) {

		int idx = 0;
		for (SAMSequenceRecord ssr : ssrList2) {

			String fname = f1.getName();
			int chridx = fname.indexOf("_");
			if (chridx < 0)
				continue;
			String chr = fname.substring(0, chridx);
			if (equals(ssr.getSequenceName(), chr)) {
				return idx;
			}
			idx++;
		}
		return idx;
	}

	
	private int comp(File f1,File f2) {

		String fname = f1.getName();
		int lastidx = fname.lastIndexOf("_");
		String ss = fname.substring(lastidx-2,lastidx);
		ss = ss.replaceAll("_","");
		
		String fname2 = f2.getName();
		int lastidx2 = fname2.lastIndexOf("_");
		String ss2 = fname2.substring(lastidx2-2,lastidx2);
		ss2 = ss2.replaceAll("_","");
		
		try{
			return Integer.parseInt(ss) - Integer.parseInt(ss2);
		}catch(Exception ex){}
		
		return ss.compareTo(ss2);
		
		
	}

	
	private boolean equals(String s1, String s2) {

		return s1.replaceAll("chr", "").equals(s2.replaceAll("chr", ""));
	}

}
