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
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.SNVHolderPlusACnv;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class GetAllelicCNVGraph {

	private final static int SIZE = 2;

	public static List<DisplayObject> getChartLists(AllelicCNV alCNV, String id) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//
		list.add(new DisplayObject(cnvList(alCNV), SIZE, "Allelic CNV "+id));
		return list;
	}

	private static JFreeChart cnvList(AllelicCNV alCNV) {

		int tsize = 0;
		int cnt = 0;
		XYSeries series1 = new XYSeries("allele1 row");
		XYSeries series2 = new XYSeries("allele2 row");

		XYSeries series3 = new XYSeries("allele1 gc correction");
		XYSeries series4 = new XYSeries("allele2 gc correction");

		XYSeries series5 = new XYSeries("allele1 wt");
		XYSeries series6 = new XYSeries("allele2 wt");

		XYSeries series7 = new XYSeries("allele1");
		XYSeries series8 = new XYSeries("allele2");

		List<Integer> chrMark = new ArrayList<Integer>();
		chrMark.add(0);

		for (List<SNVHolderPlusACnv> plist : alCNV.getList()) {
			tsize = tsize + plist.size();
			chrMark.add(tsize);
			for (SNVHolderPlusACnv sc : plist) {
				cnt++;
				series1.add(cnt, sc.getHighera().getRow());
				series2.add(cnt, sc.getLowera().getRow());

				series3.add(cnt, sc.getHighera().getGcadjusted());
				series4.add(cnt, sc.getLowera().getGcadjusted());

				series5.add(cnt, sc.getHighera().getWtval());
				series6.add(cnt, sc.getLowera().getWtval());

				series7.add(cnt, sc.getHighera().getHmmval());
				series8.add(cnt, sc.getLowera().getHmmval());
			}
		}

		// //
		// make a common vertical axis for all the sub-plots
		NumberAxis xAxis = new NumberAxis("pos");
		// xAxis.setRange(0, tsize);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis("n");
		yAxis.setRange(-1, 5);
		yAxis.setAutoRangeIncludesZero(false);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.RED);
		renderer0.setSeriesShapesVisible(1, false);
		renderer0.setSeriesPaint(1, ChartColor.BLUE);

		CombinedDomainXYPlot parent = new CombinedDomainXYPlot(xAxis);

		XYSeriesCollection data0 = new XYSeriesCollection();
		data0.addSeries(series1); // add subplot 1...
		data0.addSeries(series2); // add subplot 1...
		XYPlot subplot1 = new XYPlot(data0, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);
		addChrMaker(subplot1, chrMark);
		parent.add(subplot1, 1);

		XYSeriesCollection data1 = new XYSeriesCollection();
		data1.addSeries(series3); // add subplot 1...
		data1.addSeries(series4); // add subplot 1...
		XYPlot subplot2 = new XYPlot(data1, xAxis, yAxis, renderer0);
		subplot2.setDomainCrosshairVisible(true);
		subplot2.setRangeCrosshairVisible(true);
		addChrMaker(subplot2, chrMark);
		parent.add(subplot2, 1);
//
		XYSeriesCollection data2 = new XYSeriesCollection();
		data2.addSeries(series5); // add subplot 1...
		data2.addSeries(series6); // add subplot 1...
		XYPlot subplot3 = new XYPlot(data2, xAxis, yAxis, renderer0);
		subplot3.setDomainCrosshairVisible(true);
		subplot3.setRangeCrosshairVisible(true);
		addChrMaker(subplot3, chrMark);
		parent.add(subplot3, 1);

		XYSeriesCollection data3 = new XYSeriesCollection();
		data3.addSeries(series7); // add subplot 1...
		data3.addSeries(series8); // add subplot 1...
		XYPlot subplot4 = new XYPlot(data3, xAxis, yAxis, renderer0);
		subplot4.setDomainCrosshairVisible(true);
		subplot4.setRangeCrosshairVisible(true);
		addChrMaker(subplot4, chrMark);
		parent.add(subplot4, 1);

		// createLineChart(data, string, string2, string3);
		JFreeChart chart = new JFreeChart("cnvlist",
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;

	}

	private static void addChrMaker(XYPlot xyplot, List<Integer> chrMark) {

		for (int cnt = 0; cnt + 1 < chrMark.size(); cnt++) {
			Marker marker = new IntervalMarker(chrMark.get(cnt),
					chrMark.get(cnt + 1));
			marker.setLabel((""+cnt+1));
			marker.setLabelPaint(Color.BLACK);
			if (cnt % 2 == 0) {
				marker.setPaint(Color.LIGHT_GRAY.brighter().brighter());
			}
			marker.setAlpha(0.1f);
			xyplot.addDomainMarker(marker, Layer.BACKGROUND);
		}

	}

}
