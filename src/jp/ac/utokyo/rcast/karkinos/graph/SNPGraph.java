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

import java.awt.geom.Ellipse2D;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.jfree.chart.ChartColor;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.urls.StandardXYURLGenerator;
import org.jfree.chart.urls.XYURLGenerator;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.lowagie.text.Table;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.DepthCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;

public class SNPGraph {

	public static List<DisplayObject> getChartList(DataSet dataset) {

		//
		List<DisplayObject> list = new ArrayList<DisplayObject>();

		
		XYSeriesCollection data0 = new XYSeriesCollection();
		XYSeries series0 = null;
		Map<Double,XYSeries> md = new TreeMap<Double,XYSeries>();
		//
		for (SNVHolder snv : dataset.getSnvlist()) {

			if (snv.getDbSNPbean() != null) {

				boolean ed1 = (snv.getNormal().getTotalcnt() >= 10);
				boolean ed2 = (snv.getTumor().getTotalcnt() >= 10);
				//
				double cn = snv.getCi().getVaridateVal();				
				if(md.containsKey(cn)){
					series0 = md.get(cn);
				}else{
					series0 = new XYSeries((float)cn);
					md.put(cn, series0);					
				}	
				////
					if (ed1 && ed2) {
						DbSNPBean dbb = snv.getDbSNPbean();
						if (dbb.getMode() == DbSNPAnnotation.MODEdbSNP) {

							double x = snv.getNormal().getRatio();
							double y = snv.getTumor().getRatio();
							series0.add(x, y);								
					}					
					
				}

			}
		}	
		
		Iterator<Double> ite = md.keySet().iterator();
		while(ite.hasNext()){
			data0.addSeries(md.get(ite.next()));
		}		
		
		JFreeChart chart0 = createScatterPlot("all SNP allele frequency", "Normal",
				"Tumor", data0, PlotOrientation.VERTICAL, false, false, false);

		
		
		XYSeriesCollection data = new XYSeriesCollection();
		XYSeries series = new XYSeries("tumor/normal SNP");

		long snpRecurrent = 0;
		long snpDiff = 0;

		for (SNVHolder snv : dataset.getSnvlist()) {

			if (snv.getDbSNPbean() != null) {

				boolean ed1 = (snv.getNormal().getTotalcnt() >= 10);
				boolean ed2 = (snv.getTumor().getTotalcnt() >= 10);
				if (snv.getCi().getVaridateVal() == 1) {

					if (ed1 && ed2) {
						DbSNPBean dbb = snv.getDbSNPbean();
						if (dbb.getMode() == DbSNPAnnotation.MODEdbSNP) {

							double x = snv.getNormal().getRatio();
							double y = snv.getTumor().getRatio();
							if ((x >= 0.4)) {
								series.add(x, y);
								if (Math.abs(x - y) <= 0.35) {
									snpRecurrent++;
								} else {
									snpDiff++;
								}
							}

						}
					}

				}

			}
		}	

		data.addSeries(series);
		JFreeChart chart = createScatterPlot("SNP allele frequency", "Normal",
				"Tumor", data, PlotOrientation.VERTICAL, false, false, false);

		List oList = new ArrayList();
		oList.add(chart0);
		oList.add(chart);
		oList.add(getTable(snpRecurrent, snpDiff));

		list.add(new DisplayObject(oList, 2, "dbSNP info"));

		return list;

	}

	private static Object getTable(long snpRecurrent, long snpDiff) {
		Table table;
		try {
			table = new Table(4);
			// row1
			table.addCell("total dbSNP detected in 2n");
			table.addCell("recurrent SNP");
			table.addCell("different SNP");
			table.addCell("% diffrerence SNP");

			table.addCell(format(snpRecurrent + snpDiff));
			table.addCell(format(snpRecurrent));
			table.addCell(format(snpDiff));
			double ratio = ((double) snpDiff / (double) (snpRecurrent + snpDiff)) * 100;
			table.addCell(format(ratio)+"%");

			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static String format(long num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	public static JFreeChart createScatterPlot(String title, String xAxisLabel,
			String yAxisLabel, XYDataset dataset, PlotOrientation orientation,
			boolean legend, boolean tooltips, boolean urls) {

		if (orientation == null) {
			throw new IllegalArgumentException("Null 'orientation' argument.");
		}

		NumberAxis xAxis = new NumberAxis(xAxisLabel);
		xAxis.setRange(0, 1);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis(yAxisLabel);
		yAxis.setRange(0, 1);
		yAxis.setAutoRangeIncludesZero(false);

		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, null);

		XYToolTipGenerator toolTipGenerator = null;
		if (tooltips) {
			toolTipGenerator = new StandardXYToolTipGenerator();
		}

		XYURLGenerator urlGenerator = null;
		if (urls) {
			urlGenerator = new StandardXYURLGenerator();
		}
		XYItemRenderer renderer = new XYLineAndShapeRenderer(false, true);
		renderer.setSeriesPaint(0, ChartColor.DARK_GRAY);
		renderer.setSeriesShape(0, new Ellipse2D.Float(1.0f, 1.0f, 1.0f, 1.0f));

		renderer.setBaseToolTipGenerator(toolTipGenerator);
		renderer.setURLGenerator(urlGenerator);
		plot.setRenderer(renderer);
		plot.setOrientation(orientation);

		JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT,
				plot, legend);
		return chart;
 
	}

}
