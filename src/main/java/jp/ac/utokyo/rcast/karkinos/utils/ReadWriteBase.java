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
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;

import java.io.File;

public abstract class ReadWriteBase {

	public static SAMFileReader getReader(File INPUT) {
		
//		boolean validextension = INPUT.getName().endsWith("sam")||INPUT.getName().endsWith("bam");
//		if(!validextension){
//			
//		}
		SAMFileReader reader = new SAMFileReader(INPUT);
		SAMFileHeader sfh = reader.getFileHeader();
		if (sfh.getAttribute("SO") == null
				|| sfh.getAttribute("SO").equals("sorted")) {
			sfh.setSortOrder(SortOrder.coordinate);
		}
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		return reader;
	}

	public static SAMFileReader getReader(String in) {

		File INPUT = new File(in);
		return getReader(INPUT);
	}

}
