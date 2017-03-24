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
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.karkinos.noisefilter.AFDepthMatrix;
import jp.ac.utokyo.karkinos.noisefilter.BinData;
import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.karkinos.ploidy.MatchMatrixBean;
import jp.ac.utokyo.karkinos.ploidy.TheoreticalNodes;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.MatrixSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYZDataset;

public class GetCNVPeaksCharts {

	public static List<DisplayObject> getChartLists(PeaksInfo pi) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//

		JFreeChart jfc = getPeakPick(pi);
		list.add(new DisplayObject(jfc, 2, "CNV peaks estimates"));

		try {

			JFreeChart mattixmatch2 = getMattixMatch(pi, 2);
			list.add(new DisplayObject(mattixmatch2, 2, "matrix match"));
			JFreeChart mattixmatch4 = getMattixMatch(pi, 4);
			list.add(new DisplayObject(mattixmatch4, 2, "matrix match"));

		} catch (Exception ex) {
		}
		return list;

	}

	private static JFreeChart getMattixMatch(PeaksInfo pi, int baseploidy) {

		//
		MatchMatrixBean mmb = pi.getMatchmatrix();
		XYZDataset xyz = getMSC(mmb, baseploidy);
		int purity = mmb.getBestmme().getPurity();
		String infostr = "base ploidy=" + baseploidy;
		if(mmb.getPloidyflg() == baseploidy){
			infostr = infostr + " purity = "+ purity + " numhit=" + mmb.getBestmme().getNumhit()
			+" "+ mmb.getBestmme().getHitnodes();
		}
		JFreeChart jfc = ChartFactory.createBubbleChart(infostr,
				"AF", "peakdist", xyz, PlotOrientation.VERTICAL, true, false,
				false);
		
				
		XYPlot xyp = jfc.getXYPlot();
		xyp.setBackgroundPaint(Color.WHITE);
		
		XYItemRenderer renderer0 = xyp.getRenderer();
		renderer0.setSeriesPaint(0, ChartColor.DARK_GRAY);
		
		
		
		NumberAxis xAxis = new NumberAxis("AF");
		xAxis.setRange(0, 1);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis("pos");
		yAxis.setRange(-2, 2);
		yAxis.setAutoRangeIncludesZero(false);
		
		return jfc;

	}

	private static XYZDataset getMSC(MatchMatrixBean mmb, int baseploidy) {

		//

		List<TheoreticalNodes> tnode = null;
		if (baseploidy == 2) {
			tnode = mmb.getDiplist();
		} else {
			tnode = mmb.getTetraplist();
		}
		XYZDataset xyz = new MyXYZDataset(tnode, mmb.getList(), baseploidy);
		return xyz;

	}

	private static JFreeChart getPeakPick(PeaksInfo pi) {
		XYSeriesCollection data = new XYSeriesCollection();
		XYSeries series01 = new XYSeries("row");
		XYSeries series02 = new XYSeries("smoothed");

		XYSeriesCollection data2 = new XYSeriesCollection();
		XYSeries series05 = new XYSeries("peak analysis");

		XYSeriesCollection data1 = new XYSeriesCollection();
		XYSeries series03 = new XYSeries("EM approximation");
		XYSeries series04 = new XYSeries("additional peak");

		for (int n = 0; n < pi.getSignalcount().length; n++) {

			double x = n * 0.001;
			int rowcount = pi.getSignalcount()[n];
			double soothed = pi.getMa()[n];
			series01.add(x, rowcount);
			series02.add(x, soothed);

		}

		for (int n = 0; n < pi.getPeaksignals().length; n++) {

			double x = n * 0.001;
			double ps = pi.getPeaksignals()[n];
			series05.add(x, ps);

		}

		data.addSeries(series02);
		data.addSeries(series01);
		data1.addSeries(series03);
		data1.addSeries(series04);
		data2.addSeries(series05);

		for (double x = 0; x < 4; x = x + 0.001) {

			double y = pi.getVal(x);
			series03.add(x, pi.getAcutualVal(x));
			series04.add(x, pi.getArtifitialVal(x));

		}

		JFreeChart jfc = getGraph(data, data1, data2, PlotOrientation.VERTICAL);
		return jfc;
	}

	private static JFreeChart getGraph(XYDataset dataset0, XYDataset dataset1,
			XYSeriesCollection data2, PlotOrientation orientation) {

		NumberAxis xAxis = new NumberAxis("T/N depth ratio");
		xAxis.setRange(0, 4);
		CombinedDomainXYPlot parent = new CombinedDomainXYPlot(xAxis);

		JFreeChart chart = ChartFactory.createScatterPlot("T/N distribution",
				"T/N ratio", "count", dataset0, PlotOrientation.VERTICAL, true,
				false, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setDomainCrosshairVisible(false);
		plot.setRangeCrosshairVisible(false);
		plot.setBackgroundPaint(Color.WHITE);

		final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesShapesVisible(0, true);
		renderer1.setSeriesShapesVisible(1, false);
		// renderer1.setSeriesLinesVisible(0, false);
		plot.setRenderer(renderer1);
		parent.add(plot, 1);

		//
		JFreeChart chart3 = ChartFactory.createScatterPlot("peaksignals",
				"T/N ratio", "peaksignal", data2, PlotOrientation.VERTICAL,
				true, false, false);
		XYPlot plot3 = (XYPlot) chart3.getPlot();
		plot3.setDomainCrosshairVisible(false);
		plot3.setRangeCrosshairVisible(false);
		plot3.setBackgroundPaint(Color.WHITE);
		plot3.setDomainAxis(xAxis);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, true);
		renderer0.setSeriesLinesVisible(0, false);
		renderer0.setSeriesShapesVisible(1, true);
		renderer0.setSeriesLinesVisible(1, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setSeriesPaint(1, ChartColor.RED);
		plot3.setRenderer(0, renderer0);
		parent.add(plot3, 1);
		//
		JFreeChart chart2 = ChartFactory.createScatterPlot("T/N distribution",
				"T/N ratio", "approx func", dataset1, PlotOrientation.VERTICAL,
				true, false, false);
		XYPlot plot2 = (XYPlot) chart2.getPlot();
		plot2.setDomainCrosshairVisible(false);
		plot2.setRangeCrosshairVisible(false);
		plot2.setBackgroundPaint(Color.WHITE);
		plot2.setDomainAxis(xAxis);

		final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
		renderer2.setSeriesShapesVisible(0, true);
		renderer2.setSeriesLinesVisible(0, false);
		renderer2.setSeriesShapesVisible(1, true);
		renderer2.setSeriesLinesVisible(1, false);
		renderer2.setSeriesPaint(0, ChartColor.BLUE);
		renderer2.setSeriesPaint(1, ChartColor.RED);
		plot2.setRenderer(0, renderer2);
		parent.add(plot2, 1);

		JFreeChart chart0 = new JFreeChart("AF/depth plot",
				JFreeChart.DEFAULT_TITLE_FONT, parent, false);

		return chart0;

	}

}
