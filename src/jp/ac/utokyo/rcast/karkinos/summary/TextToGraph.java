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

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;

import com.lowagie.text.Document;
import com.lowagie.text.Image;
import com.lowagie.text.Paragraph;
import com.lowagie.text.pdf.PdfWriter;

public class TextToGraph {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// String s =
		// "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/TCGA-MR-A520-01A-11D-A25V-10.vs-10_cnvAllelicDepth.txt";
		// String s2 =
		// "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/TCGA-MR-A520-01A-11D-A25V-10.vs-10_cnvdepth.txt";
		// String sout = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/out.pdf";

		String s = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/HCC-JP-468-T.vs-N_cnvAllelicDepth.txt";
		String s2 = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/HCC-JP-468-T.vs-N_cnvdepth.txt";
		String sout = "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/468out.pdf";

		Document document = null;
		PdfWriter writer = null;
		//
		try {

			int width = 1200;
			int hight = 500;
			int size = 1;
			FileOutputStream fileOutputStream = new FileOutputStream(sout);
			document = new Document();

			writer = PdfWriter.getInstance(document, fileOutputStream);
			document.open();

			JFreeChart chart0 = getChart(new File(s2), 0);
			JFreeChart chart1 = getChart(new File(s2), 1);
			JFreeChart chart2 = getChart(new File(s2), 2);
			JFreeChart chart3 = getChart(new File(s2), 3);
			JFreeChart chart4 = getChart(new File(s), 4);
			JFreeChart chart5 = getChart(new File(s), 5);
			JFreeChart chart6 = getChart(new File(s), 6);
			
			BufferedImage bufferedImage = chart0.createBufferedImage(width
					* size, hight * size);
			Image image = Image.getInstance(writer, bufferedImage, 1.0f);
			image.scalePercent(20);
			document.add(image);

			BufferedImage bufferedImage1 = chart1.createBufferedImage(width
					* size, hight * size);
			Image image1 = Image.getInstance(writer, bufferedImage1, 1.0f);
			image1.scalePercent(20);
			document.add(image1);

			BufferedImage bufferedImage2 = chart2.createBufferedImage(width
					* size, hight * size);
			Image image2 = Image.getInstance(writer, bufferedImage2, 1.0f);
			image2.scalePercent(20);
			document.add(image2);

			BufferedImage bufferedImage3 = chart3.createBufferedImage(width
					* size, hight * size);
			Image image3 = Image.getInstance(writer, bufferedImage3, 1.0f);
			image3.scalePercent(20);
			document.add(image3);

			BufferedImage bufferedImage4 = chart4.createBufferedImage(width
					* size, hight * size);
			Image image4 = Image.getInstance(writer, bufferedImage4, 1.0f);
			image4.scalePercent(20);
			document.add(image4);

			BufferedImage bufferedImage5 = chart5.createBufferedImage(width
					* size, hight * size);
			Image image5 = Image.getInstance(writer, bufferedImage5, 1.0f);
			image5.scalePercent(20);
			document.add(image5);
			
			BufferedImage bufferedImage6 = chart6.createBufferedImage(width
					* size, hight * size);
			Image image6 = Image.getInstance(writer, bufferedImage6, 1.0f);
			image6.scalePercent(20);
			document.add(image6);

			document.close();
			document = null;
			writer.close();
			writer = null;

		} catch (Exception e) {
			e.printStackTrace();
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

	private static JFreeChart getChart(File f, int flg)
			throws NumberFormatException, IOException {

		int tsize = 0;
		int cnt = 0;
		XYSeries series1 = new XYSeries("high");
		XYSeries series2 = new XYSeries("low");

		BufferedReader br = new BufferedReader(new FileReader(f));

		String s = null;
		int count = 0;
		long pos = 0;
		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			String[] sa = s.split("\t");
			String chr0 = sa[0];
			int n = 0;
			try {
				n = Integer.parseInt(chr0);

			} catch (Exception ex) {
				break;
			}
			pos++;

			long start0 = Integer.parseInt(sa[1]);
			// int end0 = Integer.parseInt(sa[2]);
			// long pos = n * 1000000000 + start0;
			// long pos = start0;
			if (flg == 0) {
				series1.add(pos, Double.parseDouble(sa[5]));
			} else if (flg == 1) {
				series1.add(pos, Double.parseDouble(sa[6]));
			} else if (flg == 2) {
				series1.add(pos, Double.parseDouble(sa[7]));
			} else if (flg == 3) {
				series1.add(pos, Double.parseDouble(sa[8]));
				series2.add(pos, Double.parseDouble(sa[9]));
			} else if (flg == 4) {
				series1.add(pos, Double.parseDouble(sa[5]));
				series2.add(pos, Double.parseDouble(sa[6]));
			} else if (flg == 5) {

				series1.add(pos, Double.parseDouble(sa[7]));
				series2.add(pos, Double.parseDouble(sa[8]));
			} else if (flg == 6) {

				series1.add(pos, Double.parseDouble(sa[9]));
				series2.add(pos, Double.parseDouble(sa[10]));
			}
		}

		// //
		// make a common vertical axis for all the sub-plots
		NumberAxis xAxis = new NumberAxis("pos");
		NumberAxis yAxis = new NumberAxis("n");
		yAxis.setRange(0, 4);
		yAxis.setAutoRangeIncludesZero(false);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);

		XYSeriesCollection data0 = new XYSeriesCollection();
		data0.addSeries(series1); // add subplot 1...
		data0.addSeries(series2); // add subplot 1...
		XYPlot subplot1 = new XYPlot(data0, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);

		JFreeChart chart = createXYLineChart("", "", "", data0,
				PlotOrientation.VERTICAL);

		// createLineChart(data, string, string2, string3);
		// JFreeChart chart = new JFreeChart("cnvlist",
		// JFreeChart.DEFAULT_TITLE_FONT, subplot1, true);

		return chart;

	}

	public static JFreeChart createXYLineChart(String title, String xAxisLabel,
			String yAxisLabel, XYDataset dataset, PlotOrientation orientation) {

		if (orientation == null) {
			throw new IllegalArgumentException("Null 'orientation' argument.");
		}
		NumberAxis xAxis = new NumberAxis(xAxisLabel);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis(yAxisLabel);
		yAxis.setRange(0, 6);
		yAxis.setAutoRangeIncludesZero(false);
		XYItemRenderer renderer = new XYLineAndShapeRenderer(true, false);
		renderer.setSeriesPaint(0, ChartColor.BLUE);

		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
		plot.setOrientation(orientation);
		plot.setBackgroundPaint(Color.WHITE);

		JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT,
				plot, false);

		return chart;

	}

}
