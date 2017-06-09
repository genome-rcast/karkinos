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
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

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

	public SAMFileReader getReader(InputStream INPUT) {

		SAMFileReader reader = new SAMFileReader(INPUT);
		SAMFileHeader sfh = reader.getFileHeader();
		if (sfh.getAttribute("SO") == null
				|| sfh.getAttribute("SO").equals("sorted")) {
			sfh.setSortOrder(SortOrder.coordinate);
		}
		reader.setValidationStringency(ValidationStringency.SILENT);
		return reader;
	}

	public SAMFileReader getIsReader(String in) throws FileNotFoundException {

		File INPUT = new File(in);
		FileInputStream fis = new FileInputStream(INPUT);
		return getReader(fis);
	}

	public static SAMFileReader getReader(String in) {

		File INPUT = new File(in);
		return getReader(INPUT);
	}

	public static SAMFileWriter getWriter(SAMFileHeader sfh2, String OUTPUT)
			throws Exception {

		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		sfh2.setSortOrder(SortOrder.unsorted);
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(sfh2, false,
				new File(OUTPUT));

		return writer;
	}

	public static SAMFileWriter getPreSortWriter(SAMFileHeader sfh2, String OUTPUT)
			throws Exception {

		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		sfh2.setSortOrder(SortOrder.coordinate);
		factory.setCreateIndex(true);
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(sfh2, true, new File(
				OUTPUT));

		return writer;
	}

	public static SAMFileWriter getSortWriter(SAMFileHeader sfh2, String OUTPUT)
			throws Exception {

		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		sfh2.setSortOrder(SortOrder.coordinate);
		factory.setCreateIndex(true);
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(sfh2, false,
				new File(OUTPUT));

		return writer;
	}

}
