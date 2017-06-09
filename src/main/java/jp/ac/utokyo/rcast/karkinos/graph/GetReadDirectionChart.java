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
import java.util.Collection;
import java.util.List;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import jp.ac.utokyo.rcast.karkinos.bean.BaitSampling;
import jp.ac.utokyo.rcast.karkinos.distribution.XYSeriesExtract;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class GetReadDirectionChart {

	private final static int SIZE = 2;

	public static Collection<? extends DisplayObject> getChartLists(
			DataSet dataset) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		List olist = new ArrayList();
		JFreeChart chartnormal = getDirectionChart(dataset,true);
		olist.add(chartnormal);
		JFreeChart charttumor = getDirectionChart(dataset,false);
		olist.add(charttumor );
		
		list.add(new DisplayObject(olist, SIZE, "Bait sampling"));

		return list;

	}

	private static JFreeChart getDirectionChart(DataSet dataset, boolean normal) {

		BaitSampling bsnormal = dataset.getNormal().getBs();
		BaitSampling bstumor = dataset.getTumor().getBs();
		//

		NumberAxis xAxis = new NumberAxis("position in bait");
		xAxis.setRange(0, BaitSampling.bait_sample_length);
		NumberAxis yAxis = new NumberAxis("forward/reverse ratio");
		yAxis.setRange(0, 1);
		yAxis.setAutoRangeIncludesZero(true);
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);
		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, true);
		renderer0.setSeriesLinesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setSeriesPaint(1, ChartColor.RED);
		
		renderer0.setSeriesPaint(2, ChartColor.LIGHT_GRAY);
		renderer0.setSeriesPaint(3, ChartColor.LIGHT_GRAY);
		renderer0.setSeriesPaint(4, ChartColor.LIGHT_GRAY);
		renderer0.setSeriesPaint(5, ChartColor.LIGHT_GRAY);
		
		BaitSampling bs = normal ? bsnormal : bstumor;
		XYPlot subplot1 = new XYPlot(getXYSeries(bs, true), xAxis, yAxis,
				renderer0);
		XYPlot subplot2 = new XYPlot(getXYSeries(bs, false), xAxis, yAxis,
				renderer0);
		parent.add(subplot1, 1);
		parent.add(subplot2, 1);
		// createLineChart(data, string, string2, string3);
		JFreeChart chart = new JFreeChart("Bait sampling Normal",
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;
	}

	private static XYDataset getXYSeries(BaitSampling bs, boolean forward) {

		XYSeries seriesf = new XYSeries("forward rate");
		XYSeries seriesfsd1 = new XYSeries("forward rate sd+");
		XYSeries seriesfsd2 = new XYSeries("forward rate sd-");
		XYSeries seriesr = new XYSeries("reverse rate");
		XYSeries seriesrsd1 = new XYSeries("reverse rate sd+");
		XYSeries seriesrsd2 = new XYSeries("reverse rate sd-");

		for (int n = 0; n < BaitSampling.bait_sample_length; n++) {
			double fratio =0;
			double fratiosd1 =0;
			double fratiosd2 =0;
			
			double rratio = 0;
			double rratiosd1 = 0;
			double rratiosd2 = 0;
			
			if (forward) {

				 fratio = bs.getS_fstat()[n].getMean();
				 double sd = bs.getS_fstat()[n].getStandardDeviation();
				 fratiosd1 = fratio+sd;
				 fratiosd2 = fratio-sd;
				 rratio = bs.getS_rstat()[n].getMean();
				 double sd2 = bs.getS_rstat()[n].getStandardDeviation();
				 rratiosd1 = rratio+sd2;
				 rratiosd2 = rratio-sd2;
				 
			}else{
				
				 fratio = bs.getE_fstat()[n].getMean();
				 double sd = bs.getE_fstat()[n].getStandardDeviation();
				 fratiosd1 = fratio+sd;
				 fratiosd2 = fratio-sd;
				 rratio = bs.getE_rstat()[n].getMean();
				 double sd2 = bs.getE_rstat()[n].getStandardDeviation();
				 rratiosd1 = rratio+sd2;
				 rratiosd2 = rratio-sd2;
				 
			}
			
			seriesf.add(n, fratio);
			seriesr.add(n, rratio);
			
			seriesfsd1.add(n,fratiosd1);
			seriesfsd2.add(n,fratiosd2);
			seriesrsd1.add(n,rratiosd1);
			seriesrsd2.add(n,rratiosd2);
			
		}
		XYSeriesCollection data = new XYSeriesCollection();
		data.addSeries(seriesf);
		data.addSeries(seriesr);
		
		data.addSeries(seriesfsd1);
		data.addSeries(seriesfsd2);
		data.addSeries(seriesrsd1);
		data.addSeries(seriesrsd2);

		return data;
	}

}
