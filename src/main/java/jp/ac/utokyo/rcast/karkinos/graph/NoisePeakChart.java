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
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.karkinos.noisefilter.AFDepthMatrix;
import jp.ac.utokyo.karkinos.noisefilter.BinData;
import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class NoisePeakChart {
	
	public static List<DisplayObject> getChartLists(NoiseAnalysis na,double tc) {
		
		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//init filter
		double border = KarkinosProp.mintumorratio/tc;
		System.out.println("AF fixed border ="+border);
			
		AFDepthMatrix afm = na.getAfm();	
	    XYSeriesCollection data = new XYSeriesCollection();
	    XYSeriesCollection data0 = new XYSeriesCollection();
	    XYSeries series01 = new XYSeries("rejected with filter2");
	    XYSeries series02 = new XYSeries("rejected with filter1");
	    XYSeries series03 = new XYSeries("accepted");
	    
	    XYSeries series = new XYSeries("p="+KarkinosProp.pvalforNoisepeak);
		for(BinData bd:na.getBinlist()){
			
			if(bd.getDepth()<100){
				series.add(bd.getAFBorder(),bd.getDepth());
			}	
			for(Point2D p2d:bd.getAfdepth()){
				
				boolean rejected = na.reject(p2d);
				if(p2d.getX()>border){
					rejected = false;		
				}				
				if(rejected){
					series01.add(p2d.getX(),p2d.getY());
				}else{
					
					if(p2d.getX()<0.15){
						series02.add(p2d.getX(),p2d.getY());
					}else{				
						series03.add(p2d.getX(),p2d.getY());
					}	
				}
				
			}
			
		}
		data.addSeries(series);
		data0.addSeries(series03);
		data0.addSeries(series01);
		data0.addSeries(series02);
		
		
		XYSeries series2 = new XYSeries("predicted noise border");
		for(float f = 0;f<=1;f=f+0.01f){
			
			//
			series2.add(f, afm.func(f));
			
		}		
		data.addSeries(series2);
		JFreeChart jfc = getGraph(data,data0,PlotOrientation.VERTICAL);		
		list.add(new DisplayObject(jfc, 2, "noise peak limits estimates"));
		return list;
					
		
	}
	
	
	private static JFreeChart getGraph(XYDataset dataset0,XYDataset dataset1,PlotOrientation orientation) {

		NumberAxis yAxis = new NumberAxis("reads depth(TC adjusted)");
		yAxis.setRange(0, 1000);
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);
	
		
		JFreeChart chart = ChartFactory.createScatterPlot("noise peaks estimates",
				"AF (adjusted)", "depth (adjusted)",
				dataset0, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setDomainCrosshairVisible(false);
		plot.setRangeCrosshairVisible(false);
		plot.setBackgroundPaint(Color.WHITE);
		NumberAxis xAxis = new NumberAxis("AF (TC adjusted)");
		xAxis.setRange(0, 1);
		plot.setDomainAxis(xAxis);
		
		final XYLineAndShapeRenderer renderer1 = new XYLineAndShapeRenderer();
		renderer1.setSeriesShapesVisible(0, true);
		renderer1.setSeriesShapesVisible(1, false);
		//renderer1.setSeriesLinesVisible(0, false);
		plot.setRenderer(renderer1);		
		parent.add(plot, 1);
		
		JFreeChart chart2 = ChartFactory.createScatterPlot("SNV candidate",
				"AF (adjusted)", "depth (adjusted)",
				dataset1, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot2 = (XYPlot) chart2.getPlot();
		plot2.setDomainCrosshairVisible(false);
		plot2.setRangeCrosshairVisible(false);
		plot2.setBackgroundPaint(Color.WHITE);
		plot2.setDomainAxis(xAxis);
		
		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, true);
		renderer0.setSeriesLinesVisible(0, false);
		renderer0.setSeriesShapesVisible(1, true);
		renderer0.setSeriesLinesVisible(1, false);
		renderer0.setSeriesShapesVisible(2, true);
		renderer0.setSeriesLinesVisible(2, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setSeriesPaint(1, ChartColor.RED);
		plot2.setRenderer(0, renderer0);

		
		parent.add(plot2, 1);
		
		JFreeChart chart0= new JFreeChart("AF/depth plot",
				JFreeChart.DEFAULT_TITLE_FONT, parent, false);

		return chart0;


	}

}
