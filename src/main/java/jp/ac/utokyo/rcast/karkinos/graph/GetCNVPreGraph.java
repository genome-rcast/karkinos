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

import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.RegressionInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.BoxAndWhiskerToolTipGenerator;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.urls.StandardXYURLGenerator;
import org.jfree.chart.urls.XYURLGenerator;
import org.jfree.data.statistics.BoxAndWhiskerCalculator;
import org.jfree.data.statistics.BoxAndWhiskerCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class GetCNVPreGraph {

	private final static int SIZE = 2;

	public static List<DisplayObject> getChartLists(DataSet dataset) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//
		List oList = new ArrayList();
		oList.add(cnvPreList(dataset));
		//CGparcent adjust maybe implimant later
		//oList.add(cnvCGparcentBoxPlot(getBoxDataSet(dataset)));
		list.add(new DisplayObject(oList, SIZE, "CNV Analysis"));
		return list;
	}

	
	
	private static BoxAndWhiskerCategoryDataset getBoxDataSet(DataSet dataset) {
		List<List<WaveletIF>> cap = dataset.getCapInterval();
		int tsize = 0;
		int cnt = 0;
		final DefaultBoxAndWhiskerCategoryDataset ddataset 
        = new DefaultBoxAndWhiskerCategoryDataset();
		
		//devide data
		Map<Float,List<Double>> map = new TreeMap<Float,List<Double>>();
		for (List<WaveletIF> list : cap) {
			tsize = tsize + list.size();

			for (WaveletIF wi : list) {

				CapInterval civ = (CapInterval) wi;
				float cgp = civ.getCgParcent();
				double value =  wi.getValue();
				float key = round(cgp);
				List<Double> dlist = null;
				if(!map.containsKey(key)){
					dlist = new ArrayList<Double>();
					map.put(key,dlist);
				}else{
					dlist = map.get(key);
				}
				dlist.add(value);				
				cnt++;
			}
		}
		Set<Entry<Float,List<Double>>> set = map.entrySet();
		for(Entry<Float,List<Double>> entry:set){
			
			Float key = entry.getKey();
			List<Double> value = entry.getValue();
			ddataset.add(BoxAndWhiskerCalculator.calculateBoxAndWhiskerStatistics(value), key, '1');
		}		
		return ddataset;
	}



	private static float round(float originalValue) {
		 BigDecimal origin = new BigDecimal(originalValue);
        return origin.setScale(0,BigDecimal.ROUND_HALF_UP).floatValue();
        
	}



	private static JFreeChart cnvCGparcentBoxPlot(BoxAndWhiskerCategoryDataset dataset) {

	
		final CategoryAxis xAxis = new CategoryAxis("Type");
		final NumberAxis yAxis = new NumberAxis("Value");
		yAxis.setAutoRangeIncludesZero(false);
		final BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
		renderer.setFillBox(false);
		renderer.setToolTipGenerator(new BoxAndWhiskerToolTipGenerator());
		final CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis,
				renderer);

		final JFreeChart chart = new JFreeChart("CG Parcent", new Font(
				"SansSerif", Font.BOLD, 14), plot, true);

		return chart;
	}

	private static JFreeChart cnvPreList(DataSet dataset) {

		
		//prepare adjustment grapth
		XYSeries series0_0 = new XYSeries("CG Parcent Dist Median");
		XYSeries series0_1 = new XYSeries("CG Parcent regression line");
		XYSeries series0_2 = new XYSeries("CG Parcent regression line adj");
		XYSeriesCollection data1_0 = new XYSeriesCollection();
		RegressionInfo rg = dataset.getFunctionRegression().getRegInfo();
		
		for(Point2D point:rg.getGcMedianList()){
			series0_0.add(point.getX(), point.getY());
		}
		for(float f = 1;f<99;f=f+0.1f){
			
			double x = f;
			double y = rg.getReg(x);
			double y2 = rg.getRegOrg(x);
			series0_1.add(x, y);
			series0_2.add(x,y2);
		}	
		
		data1_0.addSeries(series0_0);
		data1_0.addSeries(series0_1);
		data1_0.addSeries(series0_2);
		
		
		List<List<WaveletIF>> cap = dataset.getCapInterval();
		int tsize = 0;
		int cnt = 0;
		XYSeries series1 = new XYSeries("CG Parcent");
		XYSeries series2 = new XYSeries("CG Parcent(adjusted)");
		//XYSeries series3 = new XYSeries("Capture Length");

		for (List<WaveletIF> list : cap) {
			tsize = tsize + list.size();

			for (WaveletIF wi : list) {

				CapInterval civ = (CapInterval) wi;
				series1.add(civ.getCgParcent(), civ.getOriginalValue());
				series2.add(civ.getCgParcent(), wi.getValue());
				//series2.add(civ.getLength(), wi.getValue());
				//series3.add(civ.getDuality(), civ.getOriginalValue());
				cnt++;
			}
		}
		XYSeriesCollection data1 = new XYSeriesCollection();
		data1.addSeries(series1);
		

		
		XYSeriesCollection data2 = new XYSeriesCollection();
		data2.addSeries(series2);
//		XYSeriesCollection data3 = new XYSeriesCollection();
//		data3.addSeries(series3);

		return createScatterPlotPlus(data1, data1_0, data2, 
				PlotOrientation.VERTICAL);
	}

	private static JFreeChart createScatterPlotPlus(XYDataset dataset0,
			XYDataset dataset0_1,XYDataset dataset1,  PlotOrientation orientation) {

		NumberAxis x0Axis = new NumberAxis("CG Parcent");
		NumberAxis x0_1Axis = new NumberAxis("CG Parcent Distribution Median");
		NumberAxis x1Axis = new NumberAxis("CG Parcent(adjusted)");
		//LogAxis x2Axis = new LogAxis("length");

		//NumberAxis x2Axis = new NumberAxis("capture duality");

		NumberAxis yAxis = new NumberAxis("copy number");

		yAxis.setRange(0, 5);
		yAxis.setAutoRangeIncludesZero(false);

		// make a horizontally combined plot
		// make a horizontally combined plot
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);

		XYPlot plot0 = new XYPlot(dataset0, x0Axis, yAxis, null);
		plot0.setDomainCrosshairVisible(false);
		plot0.setRangeCrosshairVisible(false);
		parent.add(plot0, 1);
		
		XYPlot plot0_1 = new XYPlot(dataset0_1, x0_1Axis, yAxis, null);
		plot0.setDomainCrosshairVisible(false);
		plot0.setRangeCrosshairVisible(false);
		parent.add(plot0_1, 1);

		XYPlot plot1 = new XYPlot(dataset1, x1Axis, yAxis, null);
		plot1.setDomainCrosshairVisible(false);
		plot1.setRangeCrosshairVisible(false);
		parent.add(plot1, 1);

//		XYPlot plot2 = new XYPlot(dataset2, x2Axis, yAxis, null);
//		plot2.setDomainCrosshairVisible(false);
//		plot2.setRangeCrosshairVisible(false);
//		parent.add(plot2, 1);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, true);
		renderer0.setSeriesLinesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		plot0.setRenderer(0, renderer0);
		plot0.setRenderer(renderer0);
		plot0.setOrientation(orientation);
		
		//dot and line
		plot0_1.setRenderer(0, renderer0);
		plot0_1.setRenderer(renderer0);
		plot0_1.setOrientation(orientation);
		final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesShapesVisible(0, true);
		renderer1.setSeriesLinesVisible(0, true);
		renderer1.setSeriesPaint(0, ChartColor.RED);
		plot0_1.setRenderer(1, renderer1);
		
		final XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
		renderer2.setSeriesShapesVisible(0, true);
		renderer2.setSeriesLinesVisible(0, true);
		renderer2.setSeriesPaint(0, ChartColor.YELLOW);
		plot0_1.setRenderer(2, renderer2);
		
		plot1.setRenderer(0, renderer0);
		plot1.setRenderer(renderer0);
		plot1.setOrientation(orientation);

//		plot2.setRenderer(0, renderer0);
//		plot2.setRenderer(renderer0);
//		plot2.setOrientation(orientation);

		JFreeChart chart = new JFreeChart("cnvanalysis",
				JFreeChart.DEFAULT_TITLE_FONT, parent, false);

		return chart;

	}

}
