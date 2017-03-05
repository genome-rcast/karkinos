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
package jp.ac.utokyo.rcast.karkinos.graph.output;

import java.awt.image.BufferedImage;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.jfree.chart.JFreeChart;

import com.lowagie.text.Annotation;
import com.lowagie.text.BadElementException;
import com.lowagie.text.Chunk;
import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Element;
import com.lowagie.text.Image;
import com.lowagie.text.Paragraph;
import com.lowagie.text.Table;
import com.lowagie.text.pdf.PdfWriter;

import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.graph.DisplayObject;
import jp.ac.utokyo.rcast.karkinos.graph.GetAllelicCNVGraph;
import jp.ac.utokyo.rcast.karkinos.graph.GetCNVCharts;
import jp.ac.utokyo.rcast.karkinos.graph.GetCNVPeaksCharts;
import jp.ac.utokyo.rcast.karkinos.graph.GetCNVPreGraph;
import jp.ac.utokyo.rcast.karkinos.graph.GetGCAdjustGrapth;
import jp.ac.utokyo.rcast.karkinos.graph.GetReadDirectionChart;
import jp.ac.utokyo.rcast.karkinos.graph.GetReadsStatsChart;
import jp.ac.utokyo.rcast.karkinos.graph.GetSNVChart;
import jp.ac.utokyo.rcast.karkinos.graph.GetSNVDistChart;
import jp.ac.utokyo.rcast.karkinos.graph.GetThresholdsCharts;
import jp.ac.utokyo.rcast.karkinos.graph.NoisePeakChart;
import jp.ac.utokyo.rcast.karkinos.graph.SNPGraph;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;

public class PdfReport {

	public static void report(ReadsSummary readsSummary, String readsStat,
			DataSet dataset, AllelicCNV alCNV, NoiseAnalysis na, 
			PeaksInfo pi,String id, String outfile) throws Exception {

		List<DisplayObject> chartList = new ArrayList<DisplayObject>();
		chartList.addAll(GetReadsStatsChart.getChartLists(readsSummary,
				readsStat));
		// chartList.addAll(GetReadDirectionChart.getChartLists(dataset));
		chartList.addAll(SNPGraph.getChartList(dataset));
		
		chartList.addAll(GetGCAdjustGrapth.getChartLists(dataset));
		chartList.addAll(GetCNVPreGraph.getChartLists(dataset));
		chartList.addAll(GetCNVPeaksCharts.getChartLists(pi));
		chartList.addAll(GetCNVCharts.getChartLists(dataset,pi,id));
		
		//allelic CNV
		chartList.addAll(GetAllelicCNVGraph.getChartLists(alCNV,id));		
		chartList.addAll(GetCNVCharts.getTables(dataset));

		// chartList.addAll(GetSNVDistChart.getChartLists(dataset));
		
		//
		try{
		 chartList.addAll(NoisePeakChart .getChartLists(na,dataset.getTumorRatio()));
		}catch(Exception ex){
			ex.printStackTrace();
		}
		chartList.addAll(GetSNVChart.getChartLists(dataset,pi));
		chartList.addAll(GetThresholdsCharts.getChartLists(dataset));

		Document document = null;
		PdfWriter writer = null;
		FileOutputStream fileOutputStream = new FileOutputStream(outfile);

		try {
			// instantiate document and writer
			document = new Document();
			writer = PdfWriter.getInstance(document, fileOutputStream);
			// open document
			document.open();
			// add image
			int width = 1200;
			int hight = 1000;

			int i = 0;
			int figcount = 0;
			for (DisplayObject dobj : chartList) {

				int size = dobj.getSize();
				Object obj = dobj.getObject();
				if (obj == null)
					continue;

				document.add(new Paragraph(String.valueOf(dobj.getTitle())));
				if (obj instanceof List) {

					for (Object childObj : (List) obj) {
						addObj(document, childObj, width, hight, size, writer,figcount);
						if ((childObj instanceof Table)
								|| (obj instanceof JFreeChart)) {
							figcount++;
							document.add(Chunk.NEXTPAGE);
						}
					}
				} else {
					addObj(document, obj, width, hight, size, writer,figcount);
					figcount++;
				}
				document.add(Chunk.NEXTPAGE);
				i++;
			}
			// release resources
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

	private static void addObj(Document document, Object obj, int width,
			int hight, int size, PdfWriter writer, int figcount) throws Exception,
			IOException {
		if (obj instanceof JFreeChart) {

			JFreeChart chart = (JFreeChart) obj;
			BufferedImage bufferedImage = chart.createBufferedImage(width
					* size, hight * size);
			Image image = Image.getInstance(writer, bufferedImage, 1.0f);
			image.scalePercent(20);
			image.setAnnotation(new Annotation("fig"+figcount, "karkinos fig"+figcount));
			document.add(image);

		} else if (obj instanceof Element) {

			document.add((Element) obj);

		} else if (obj instanceof java.awt.Image) {

			java.awt.Image aim = (java.awt.Image)obj;
			Image image = Image.getInstance(writer, aim, 1.0f);
			image.scalePercent(50);
			document.add(image);

		} else {

			document.add(new Paragraph(String.valueOf(obj)));

		}

	}

}
