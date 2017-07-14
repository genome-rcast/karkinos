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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.IntervalXYDataset;

import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.distribution.DataHolderByCN;
import jp.ac.utokyo.rcast.karkinos.distribution.XYSeriesExtract;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class GetSNVDistChart {
	
	private final static int SIZE = 1;
	
	public static List<DisplayObject> getChartLists(DataSet dataset) {

		AnalyseDist analyseDist = dataset.getAnalyseDist();
		Map<Float, DataHolderByCN> snvMAP = analyseDist.getMap();
		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//
		XYSeriesExtract xyse = new XYSeriesExtract(snvMAP);
		List oList = new ArrayList();
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.5f) {
			IntervalXYDataset normal = xyse.getXYSeriesForLogdist(degree,true);
			IntervalXYDataset tumor = xyse.getXYSeriesForLogdist(degree, false);

			oList.add(getHetroSNPGraph(normal, tumor, degree));
		}
		list.add(new DisplayObject(oList, SIZE, "mutation Graph"));

		
		//
		return list;
	
	}	

	private static JFreeChart getHetroSNPGraph(IntervalXYDataset normal,
			IntervalXYDataset tumor, float degree) {

		int n = (int) (degree * 2);
		NumberAxis xAxis = new NumberAxis("log odds");
		xAxis.setRange(0, 100);
		//LogAxis xAxis = new LogAxis("log odds");
		//xAxis.setRange(.001, 1000);
		
		NumberAxis yAxis = new NumberAxis("counts");
		yAxis.setAutoRangeIncludesZero(true);

		// make a horizontally combined plot
		// make a horizontally combined plot
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);

		final XYItemRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesPaint(0, ChartColor.BLUE);

		// add subplot 1...

		XYPlot subplot1 = new XYPlot(normal, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);
		parent.add(subplot1, 1);

		// add subplot 2...
		final XYItemRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesPaint(0, ChartColor.BLUE);

		XYPlot subplot2 = new XYPlot(tumor, xAxis, yAxis, renderer1);
		subplot2.setDomainCrosshairVisible(true);
		subplot2.setRangeCrosshairVisible(true);
		parent.add(subplot2, 1);

		// createLineChart(data, string, string2, string3);
		JFreeChart chart = new JFreeChart("Mutation allele frequency N=" + n,
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;
	}
	
	
	
}
