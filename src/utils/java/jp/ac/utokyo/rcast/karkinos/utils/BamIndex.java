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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;



public class BamIndex extends ReadWriteBase{
	
	public static void main(String[] arg){
		
		//
		String bam ="/GLUSTER_DIST/data/users/ueda/pgtest/testbam3.bam";
		
		//

			SAMFileReader bamr = getReader(bam);
			CloseableIterator<SAMRecord> ite = bamr.iterator();
			
			SAMFileHeader sfh = bamr.getFileHeader();		
			SAMSequenceDictionary sequenceDictionary = sfh.getSequenceDictionary();
			SAMSequenceDictionary sequenceDictionarynew = new SAMSequenceDictionary();
			
			int b4 = 0;
			SAMRecord samb4 = null; 

			while(ite.hasNext()){
				
				SAMRecord sam = null;
				
				try{
					sam = ite.next();
				}catch(Exception ex){
					ex.printStackTrace();
					continue;
				}
				if(sam.getAttribute("YB")==null){
					System.out.println(sam.format());
				}
				
				//System.out.println(sam.getAttribute("YB")+" " + sam.getAttribute("YN") +" " + sam.format());		
				if(b4>sam.getAlignmentStart()){
					
					System.out.println(sam.getAttribute("YB")+" " + sam.getAttribute("YN") +" " + samb4.format());
					System.out.println(sam.getAttribute("YB")+" " + sam.getAttribute("YN") +" " + sam.format());	
				}
				
				b4 = sam.getAlignmentStart();
				samb4 = sam;
				
			}
			bamr.close();
			
		
		
	}
	
	
}
