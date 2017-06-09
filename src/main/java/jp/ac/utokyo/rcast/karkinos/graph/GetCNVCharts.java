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
package jp.ac.utokyo.rcast.karkinos.graph;

import java.awt.Color;
import java.awt.Font;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleInsets;

import com.lowagie.text.Table;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class GetCNVCharts {

	private final static int SIZE = 2;

	public static List<DisplayObject> getChartLists(DataSet dataset,
			PeaksInfo pi,String id) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//
		// list.add(cnvList(dataset));
		// list.add(hetroSNP(dataset));
		// list.add(mutation(dataset));
		list.add(new DisplayObject(cnvList(dataset, pi), SIZE, "CNV "+ id));
		return list;
	}

	public static List<DisplayObject> getTables(DataSet dataset) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		list.add(new DisplayObject(cnvSummaryTable(dataset, 1), SIZE,
				"CNV summary table"));
		list.add(new DisplayObject(cnvSummaryTable(dataset, 2), SIZE,
				"allelic CNV summary table"));
		list.add(new DisplayObject(cnvSummaryTable(dataset, 3), SIZE,
				"homozygous delation & amplification"));
		return list;
	}

	private static Object cnvSummaryTable(DataSet dataset, int dispFlg) {

		Table table;
		try {

			if (dispFlg == 1) {
				table = new Table(8);
			} else {
				table = new Table(6);
			}
			// header
			table.addCell("chr");
			table.addCell("start");
			table.addCell("end");
			table.addCell("copy number");
			table.addCell("gain/loss");
			table.addCell("num hetro SNP");

			if (dispFlg == 1) {
				table.addCell("AAF");
				table.addCell("BAF");
			}
			// table.addCell("SNP ratio correl");
			// table.addCell("rejected");

			for (CopyNumberInterval cni : dataset.getCopyNumberIntervalList(9)) {

				if (cni.getCopynumber() == dataset.getBaseploidy()){
					if(cni.getAaf()==dataset.getBaseploidy()/2){
						continue;
					}
				}
				if (dispFlg == 1) {
					if (cni.isAllelic())
						continue;
//					if (cni.getCopynumber() == 0)
//						continue;
//					if (cni.getCopynumber() >= 4)
//						continue;
				} else if (dispFlg == 2) {
					if (!cni.isAllelic())
						continue;
					if (cni.getCopynumber() == 0)
						continue;
					if (cni.isHdelation())
						continue;
				} else {
					if (1 <= cni.getCopynumber() && !cni.isHdelation()) {
						continue;
					}
					if(cni.getCopynumber()!=0 && cni.isHdelation()){
						continue;
					}
					
				}
				table.addCell(cni.getChr());
				table.addCell(format(cni.getStart()) + "");
				table.addCell(format(cni.getEnd()) + "");
				table.addCell("n=" + cni.getCopynumber());
				String s = cni.getCopynumber() < 2 ? "loss" : "gain";
				table.addCell(s);
				table.addCell(format(cni.getNoSNP()) + "");

				if (dispFlg == 1) {
					table.addCell(cni.getAaf() + "");
					table.addCell(cni.getBaf() + "");
				}

				// table.addCell(format(cni.getSnpclrrel())+"");
				// table.addCell(cni.getRejected()+"");

			}
			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	private static JFreeChart cnvList(DataSet dataset, PeaksInfo pi) {

		List<List<WaveletIF>> cap = dataset.getCapInterval();
		int tsize = 0;
		int cnt = 0;
		XYSeries series1 = new XYSeries("row");
		XYSeries series2 = new XYSeries("row(gc adjust)");
		XYSeries series3 = new XYSeries("denoise");
		// XYSeries series2_0 = new XYSeries("baselevel_1n");
		// XYSeries series2_1 = new XYSeries("baselevel_2n");
		// XYSeries series2_2 = new XYSeries("baselevel_3n");
		// XYSeries series2_3 = new XYSeries("baselevel_3n");
		// XYSeries series3 = new XYSeries("step func");
		XYSeries series4 = new XYSeries("hmm");
		XYSeries series5 = new XYSeries("colvaridated");

		XYSeries series6 = new XYSeries("AAF allelic");
		XYSeries series7 = new XYSeries("BAF allelic");

		List<Integer> chrMark = new ArrayList<Integer>();
		chrMark.add(0);
		for (List<WaveletIF> list : cap) {
			tsize = tsize + list.size();
			chrMark.add(tsize);
			for (WaveletIF wi : list) {

				if (cnt % 2 == 0) {
					// row data too heavy
					// draw half
					series1.add(cnt, ((CapInterval) wi).getOriginalValue());
				}
				if (cnt % 2 == 0) {
					// row data too heavy
					// draw half
					series2.add(cnt, wi.getValue());
				}
				series3.add(cnt, wi.getDenioseValue());
				// if(cnt%2==0){
				// series2_0.add(cnt,dataset.getBaselineLOHEstimate());
				// series2_1.add(cnt,1);
				// series2_2.add(cnt,2-dataset.getBaselineLOHEstimate());
				// series2_3.add(cnt,3-(2*dataset.getBaselineLOHEstimate()));
				// }
				// series3.add(cnt,wi.getCN());
				series4.add(cnt, wi.getHMMValue()*2);
				series5.add(cnt, ((CapInterval) wi).getVaridateVal()*2);

				series6.add(cnt, ((CapInterval) wi).getAafreq());
				series7.add(cnt, ((CapInterval) wi).getBafreq());

				cnt++;
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

		CombinedDomainXYPlot parent = new CombinedDomainXYPlot(xAxis);
		
		XYSeriesCollection data0 = new XYSeriesCollection();
		data0.addSeries(series1); // add subplot 1...
		XYPlot subplot1 = new XYPlot(data0, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);
		addChrMaker(subplot1, chrMark, dataset.getChromList());
		parent.add(subplot1, 1);

		// add subplot 2...
		final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesShapesVisible(0, false);
		renderer1.setSeriesShapesVisible(1, false);
		renderer1.setSeriesShapesVisible(2, false);
		renderer1.setSeriesShapesVisible(3, false);

		renderer1.setSeriesPaint(0, ChartColor.BLUE);
		// renderer1.setSeriesPaint(1, ChartColor.VERY_LIGHT_RED);
		// renderer1.setSeriesPaint(2, ChartColor.VERY_LIGHT_RED);
		// renderer1.setSeriesPaint(3, ChartColor.VERY_LIGHT_RED);
		// renderer1.setSeriesPaint(4, ChartColor.VERY_LIGHT_RED);

		XYSeriesCollection data1 = new XYSeriesCollection();
		data1.addSeries(series2);
		// data1.addSeries(series2_0);
		// data1.addSeries(series2_1);
		// data1.addSeries(series2_2);
		// data1.addSeries(series2_3);
		XYPlot subplot2 = new XYPlot(data1, xAxis, yAxis, renderer1);
		subplot2.setDomainCrosshairVisible(true);
		subplot2.setRangeCrosshairVisible(true);
		addMaker(subplot2, pi.getPeaklist());
		// addMaker(subplot2,dataset.getBaselineLOHEstimate());
		addChrMaker(subplot2, chrMark, dataset.getChromList());
		parent.add(subplot2, 1);

		// add subplot 3...
		final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
		renderer2.setSeriesShapesVisible(0, false);
		renderer2.setSeriesPaint(0, ChartColor.BLUE);

		XYSeriesCollection data2 = new XYSeriesCollection();
		data2.addSeries(series3);
		XYPlot subplot3 = new XYPlot(data2, xAxis, yAxis, renderer2);
		subplot3.setDomainCrosshairVisible(true);
		subplot3.setRangeCrosshairVisible(true);
		addMaker(subplot3, pi.getPeaklist());
		addChrMaker(subplot3, chrMark, dataset.getChromList());
		parent.add(subplot3, 1);

		// add subplot 4...
		final XYLineAndShapeRenderer renderer3 = new XYLineAndShapeRenderer();
		renderer3.setSeriesShapesVisible(0, false);
		renderer3.setSeriesPaint(0, ChartColor.BLUE);

		NumberAxis yAxis2 = new NumberAxis("n");
		yAxis2.setRange(0, 12);
		yAxis2.setAutoRangeIncludesZero(false);
		
		XYSeriesCollection data3 = new XYSeriesCollection();
		data3.addSeries(series4);
		XYPlot subplot4 = new XYPlot(data3, xAxis, yAxis2, renderer3);
		subplot4.setDomainCrosshairVisible(true);
		subplot4.setRangeCrosshairVisible(true);
		addChrMaker(subplot4, chrMark, dataset.getChromList());
		parent.add(subplot4, 1);

		final XYLineAndShapeRenderer renderer4 = new XYLineAndShapeRenderer();
		renderer4.setSeriesShapesVisible(0, false);
		renderer4.setSeriesPaint(0, ChartColor.BLUE);
		// subplot 4
		XYSeriesCollection data4 = new XYSeriesCollection();
		data4.addSeries(series5);
		XYPlot subplot5 = new XYPlot(data4, xAxis, yAxis2, renderer4);
		subplot5.setDomainCrosshairVisible(true);
		subplot5.setRangeCrosshairVisible(true);
		addChrMaker(subplot5, chrMark, dataset.getChromList());
		parent.add(subplot5, 1);


		final XYLineAndShapeRenderer renderer5 = new XYLineAndShapeRenderer();
		renderer5.setSeriesShapesVisible(0, false);
		renderer5.setSeriesPaint(0, ChartColor.RED);
		renderer5.setSeriesPaint(1, ChartColor.BLUE);
		// subplot 4
		XYSeriesCollection data5 = new XYSeriesCollection();
		data5.addSeries(series6);
		data5.addSeries(series7);
		XYPlot subplot6 = new XYPlot(data5, xAxis, yAxis2, renderer5);
		subplot6.setDomainCrosshairVisible(true);
		subplot6.setRangeCrosshairVisible(true);
		addChrMaker(subplot6, chrMark, dataset.getChromList());
		parent.add(subplot6, 1);

		// createLineChart(data, string, string2, string3);
		JFreeChart chart = new JFreeChart("cnvlist",
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;

	}

	private static void addChrMaker(XYPlot xyplot, List<Integer> chrMark,
			List<String> labels) {

		for (int cnt = 0; cnt + 1 < chrMark.size(); cnt++) {
			Marker marker = new IntervalMarker(chrMark.get(cnt),
					chrMark.get(cnt + 1));
			String label = "";
			if (labels.size() > cnt) {
				label = labels.get(cnt);
			}
			marker.setLabel(label);
			marker.setLabelPaint(Color.BLACK);
			if (cnt % 2 == 0) {
				marker.setPaint(Color.LIGHT_GRAY.brighter().brighter());
			}
			marker.setAlpha(0.1f);
			xyplot.addDomainMarker(marker, Layer.BACKGROUND);
		}

	}

	private static void addMaker(XYPlot xyplot, List<Peak> list) {

		for (Peak peak : list) {
			Marker marker0 = new ValueMarker(peak.getU());
			marker0.setPaint(Color.RED);
			xyplot.addRangeMarker(marker0);
		}

	}

	private static void addMaker(XYPlot xyplot, double baselineLOHEstimate) {
		Marker marker0 = new ValueMarker(baselineLOHEstimate);
		marker0.setPaint(Color.RED);
		Marker marker1 = new ValueMarker(1);
		marker1.setPaint(Color.RED);
		Marker marker2 = new ValueMarker(2 - baselineLOHEstimate);
		marker2.setPaint(Color.RED);
		Marker marker3 = new ValueMarker(3 - (2 * baselineLOHEstimate));
		marker3.setPaint(Color.RED);
		xyplot.addRangeMarker(marker0);
		xyplot.addRangeMarker(marker1);
		xyplot.addRangeMarker(marker2);
		xyplot.addRangeMarker(marker3);

	}

	//
	//
	//
	// private static DisplayObject caprate(DataSet dataset) {
	// DisplayObject dobj = new DisplayObject();
	// dobj.setObject(getChart(1,dataset));
	// dobj.setSize(SIZE);
	// return dobj;
	// }
	//
	// private static DisplayObject cnvList(DataSet dataset) {
	//
	// DisplayObject dobj = new DisplayObject();
	// dobj.setObject(getcnvlist(dataset));
	// dobj.setSize(SIZE);
	// return dobj;
	// }
	//
	// private static DisplayObject hetroSNP(DataSet dataset) {
	//
	// DisplayObject dobj = new DisplayObject();
	// dobj.setObject(getChart(2,dataset));
	// dobj.setSize(SIZE);
	// return dobj;
	// }
	//
	// private static DisplayObject mutation(DataSet dataset) {
	//
	// DisplayObject dobj = new DisplayObject();
	// dobj.setObject(getChart(3,dataset));
	// dobj.setSize(SIZE);
	// return dobj;
	//
	// }

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static String format(int num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

}
