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
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
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

import com.lowagie.text.Annotation;
import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Image;
import com.lowagie.text.Paragraph;
import com.lowagie.text.pdf.PdfWriter;

public class RegionGraph {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String indir = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_NCC/karkinosver4.0.25,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor/HCC-JP-453-T-HCC-JP-453-N/HCC-JP-453-T-HCC-JP-453-N,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/tcga,"
			+ "/GLUSTER_DIST/data/users/ueda/ICGC/bcmtest/goss,";

		//
		// String indir =
		// "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast/THCC_167T-THCC_167N/THCC_167T-THCC_167N/";
		String out = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/chr5Graprh.pdf";
		String outpng = "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/chr1Graprh.jpeg";

		List<File> list = new ArrayList<File>();

		for (String s : indir.split(",")) {
			File f = new File(s);
			// chage to reccursive file search
			searchRecursive(f, list);
		}

		 String chr = "chr5";
		 int start = 0;
		 int end = 50000000;

		// String chr = "chr11";
		// int start = 0;
		// int end = 135000000;

//		String chr = "chr1";
//		int start = 0;
//		int end = 247000000;
		int[][] target = new int[][] { { 1278756,1295162 } };

		Document document = null;
		PdfWriter writer = null;
		int width = 1200;
		int hight = 500;
		int size = 1;
		try {

			FileOutputStream fileOutputStream = new FileOutputStream(out);
			document = new Document();

			writer = PdfWriter.getInstance(document, fileOutputStream);
			document.open();

			for (File f : list) {

				document.add(new Paragraph(String.valueOf(f.getName())));
				Object[] obj = getChart(f, chr, start, end, target);
				document.add(new Paragraph("mean=" + toStr(obj[1])
						+ "\t mean original=" + toStr(obj[2])));
				JFreeChart chart = (JFreeChart) obj[0];

				// try {
				// File jpn = new File(outpng);
				// ChartUtilities.saveChartAsJPEG(jpn, chart, 300, 300);
				//
				// } catch (Exception ex) {
				//
				// }

				BufferedImage bufferedImage = chart.createBufferedImage(width
						* size, hight * size);
				Image image = Image.getInstance(writer, bufferedImage, 1.0f);
				image.scalePercent(20);
				document.add(image);

			}

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

	private static String toStr(Object obj) {

		SummaryStatistics ss0 = (SummaryStatistics) obj;
		return "n=" + ss0.getN() + ",mean =" + ss0.getMean();
	}

	private static boolean searchRecursive(File f, final List<File> list) {

		File[] listFiles = f.listFiles(new FilenameFilter() {

			public boolean accept(File dir, String name) {

				if (name.startsWith("."))
					return false;
				boolean target = dir.isDirectory()
						|| name.contains("cnvdepth.txt");
				if (!target)
					return false;
				//
				String absPath = dir.getAbsolutePath() + File.separator + name;
				if (new File(absPath).isFile()) {
					return true;
				} else {
					return searchRecursive(new File(absPath), list);
				}
			}

		});

		if (listFiles == null)
			return false;
		for (File ff : listFiles) {
			if (ff.isFile() && ff.getName().contains("cnvdepth.txt")) {

				System.out.println(ff.getName());
				list.add(ff);
			}
		}

		return true;

	}

	private static Object[] getChart(File f, String chr, int start, int end,
			int[][] target) throws NumberFormatException, IOException {

		int tsize = 0;
		int cnt = 0;
		XYSeries series1 = new XYSeries("row");
		BufferedReader br = new BufferedReader(new FileReader(f));

		SummaryStatistics ss0 = new SummaryStatistics();
		SummaryStatistics ss1 = new SummaryStatistics();

		String s = null;
		int count = 0;
		while ((s = br.readLine()) != null) {

			if (s.startsWith("#"))
				continue;
			String[] sa = s.split("\t");
			String chr0 = sa[0];
			int start0 = Integer.parseInt(sa[1]);
			int end0 = Integer.parseInt(sa[2]);
			if(!chr0.contains("chr")){
				chr0 = "chr"+chr0;
			}
			
			if (chr.equals(chr0)) {

				count++;
				series1.add(start0, Double.parseDouble(sa[5]));

				boolean intarget = interget(start0, end0, target);
				if (intarget) {
					//
					ss0.addValue(Double.parseDouble(sa[5]));
					ss1.addValue(Double.parseDouble(sa[6]));
				}

			}

		}

		// //
		// make a common vertical axis for all the sub-plots
		NumberAxis xAxis = new NumberAxis("pos");
		xAxis.setRange(0, tsize);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis("n");
		yAxis.setRange(0, 4);
		yAxis.setAutoRangeIncludesZero(false);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);

		XYSeriesCollection data0 = new XYSeriesCollection();
		data0.addSeries(series1); // add subplot 1...
		XYPlot subplot1 = new XYPlot(data0, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);

		JFreeChart chart = createXYLineChart("", "", "", data0,
				PlotOrientation.VERTICAL, target);

		// createLineChart(data, string, string2, string3);
		// JFreeChart chart = new JFreeChart("cnvlist",
		// JFreeChart.DEFAULT_TITLE_FONT, subplot1, true);

		return new Object[] { chart, ss0, ss1 };

	}

	private static boolean interget(int start0, int end, int[][] target) {
		
		for(int[] tg:target){
			
			if(tg[0]<end && tg[1] >start0){
				return true;
			}
			
		}		
		return false;
	}

	public static JFreeChart createXYLineChart(String title, String xAxisLabel,
			String yAxisLabel, XYDataset dataset, PlotOrientation orientation,
			int[][] target) {

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

		// TERT
		for (int[] tg : target) {
			Marker marker = new IntervalMarker(tg[0], tg[1]);
			marker.setOutlinePaint(Color.RED);
			plot.addDomainMarker(marker, Layer.BACKGROUND);
		}

		JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT,
				plot, false);

		return chart;

	}
}
