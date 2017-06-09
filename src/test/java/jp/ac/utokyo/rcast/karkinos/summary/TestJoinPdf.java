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
package jp.ac.utokyo.rcast.karkinos.summary;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import com.lowagie.text.Annotation;
import com.lowagie.text.Chunk;
import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Image;
import com.lowagie.text.Paragraph;
import com.lowagie.text.pdf.PRStream;
import com.lowagie.text.pdf.PdfDocument;
import com.lowagie.text.pdf.PdfName;
import com.lowagie.text.pdf.PdfObject;
import com.lowagie.text.pdf.PdfReader;
import com.lowagie.text.pdf.PdfStream;
import com.lowagie.text.pdf.PdfWriter;

public class TestJoinPdf {

	/**
	 * @param args
	 * @throws IOException
	 * @throws DocumentException
	 */
	public static void main(String[] args) throws IOException,
			DocumentException {

		File testf = new File(
				"/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/HX5T/HX5T_pdfdata.pdf");

		File testout = new File(
				"/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/HX5T/HX5T_pdfdata_exst.pdf");

		Document document = null;
		PdfWriter writer = null;
		try {

			FileOutputStream fileOutputStream = new FileOutputStream(testout);
			document = new Document();

			writer = PdfWriter
					.getInstance(document, fileOutputStream);
			document.open();
			document.add(new Paragraph("start "));
			PdfReader reader = new PdfReader(testf.getAbsolutePath());
			int xrefno = reader.getXrefSize();
			for (int n = 1; n <= xrefno; n++) {
				//
				PdfObject pdfobj = reader.getPdfObject(n);
				if (pdfobj != null) {

					if (pdfobj.isStream()) {

						PdfStream stream = (PdfStream) pdfobj;
						PdfObject pdfsubtype = stream.get(PdfName.SUBTYPE);
						if (pdfsubtype != null) {
							if (pdfsubtype.toString().equals(
									PdfName.IMAGE.toString())) {

								// image found
								PRStream prs = (PRStream)pdfobj;
								Image image = Image.getInstance(prs.getBytes());
								image.scalePercent(20);								
								document.add(image);
								
							}
						}

					}

				}
			}

			document.close();
			document = null;
			writer.close();
			writer = null;
		} catch (DocumentException de) {
			throw de;
		} catch (IOException ioe) {
			throw ioe;
		} finally {
			// release resources
			if (null != document) {
				try {
					document.close();
				} catch (Exception ex) {
				}
			}
			if (null != writer) {
				try {
					writer.close();
				} catch (Exception ex) {
				}
			}
		}

	}
}
