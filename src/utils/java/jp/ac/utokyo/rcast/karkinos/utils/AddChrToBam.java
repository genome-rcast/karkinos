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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;

public class AddChrToBam extends ReadWriteBase{
	
	
	public static void main(String[] arg) throws Exception{
		
		String bam = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/"
			+"synthetic.challenge.set4.normal.bam";
		
		String bamo = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/"
			+"synthetic.challenge.set4.normal_m.bam";
		
		SAMFileReader bamr = getReader(bam);
		CloseableIterator<SAMRecord> ite = bamr.iterator();
		
		SAMFileHeader sfh = bamr.getFileHeader();		
		SAMSequenceDictionary sequenceDictionary = sfh.getSequenceDictionary();
		SAMSequenceDictionary sequenceDictionarynew = new SAMSequenceDictionary();
		for(SAMSequenceRecord ssr:sequenceDictionary.getSequences()){
			
			String chrname = addChr(ssr.getSequenceName());
			System.out.println(chrname);
			SAMSequenceRecord sequenceRecordn = new SAMSequenceRecord(chrname,ssr.getSequenceLength());
			sequenceDictionarynew.addSequence(sequenceRecordn);
			
		}
		sfh.setSequenceDictionary(sequenceDictionarynew);
		
		
		SAMFileWriter sr = getPreSortWriter(sfh,bamo);		
		
		while(ite.hasNext()){
			
			SAMRecord sam = ite.next();
			sr.addAlignment(sam);				
			
		}
		bamr.close();
		sr.close();
		
		
		
	}

	private static String addChr(String s) {
		
		if(s.equals("X")){
			return "chrX";
		}
		if(s.equals("Y")){
			return "chrY";
		}
		if(s.equals("MT")){
			return "chrMT";
		}
		
		
		try{
		 Integer.parseInt(s);
		}catch(Exception ex){
			return s;
		}
		return "chr"+s;
	}

}
