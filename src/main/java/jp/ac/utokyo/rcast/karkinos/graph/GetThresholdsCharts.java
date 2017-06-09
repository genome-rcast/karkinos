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

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.IntervalXYDataset;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;

public class GetThresholdsCharts {

	public static Collection<? extends DisplayObject> getChartLists(
			DataSet dataset) {

		List<DisplayObject> list = new ArrayList<DisplayObject>();
		list.add(new DisplayObject(getThresHoldsChart(dataset), 1, "Thresholds"));

		return list;

	}

	private static List<JFreeChart> getThresHoldsChart(DataSet dataset) {

		List<JFreeChart> olist = new ArrayList<JFreeChart>();
		olist.add(createHistogram("adjusted ratio >"
				+KarkinosProp.mintumorratio,
				 "ratio", "count", dataset,
				1,KarkinosProp.mintumorratio));

		olist.add(createHistogram("original ratio >"
				+ KarkinosProp.min_initial_tumorratio, "ratio", "count", dataset, 2,
				KarkinosProp.min_initial_tumorratio));

		olist.add(createHistogram("rms of base Quality", "BQ", "count",
				dataset, 3));

		olist.add(createHistogram("rms of mapping Quality>"
				+ KarkinosProp.minMapQ, "MQ", "count", dataset, 4,
				KarkinosProp.minMapQ));

		olist.add(createHistogram("sum phred quality >"
				+ KarkinosProp.minPhredQual, "Phred Qual", "count", dataset, 5,
				KarkinosProp.minPhredQual));

		olist.add(createHistogram("Seq entropy >" + KarkinosProp.minEntropy,
				"Seq entropy", "count", dataset, 6, KarkinosProp.minEntropy));

		olist.add(createHistogram(
				"Mappability >" + KarkinosProp.minMappability, "mappability",
				"count", dataset, 7, KarkinosProp.minMappability));

		olist.add(createHistogram("pval for directional reads >"
				+ KarkinosProp.Fisher_Thres_For_Reads_Direction,
				"pval for direc", "count", dataset, 8,
				KarkinosProp.Fisher_Thres_For_Reads_Direction));

		olist.add(createHistogram("pval tumor snv fisher <"
				+ KarkinosProp.Fisher_Thres_For_SNV_Detection, "pval", "count",
				dataset, 9, KarkinosProp.Fisher_Thres_For_SNV_Detection));

		olist.add(createHistogram("log odds normal", "odds score", "count",
				dataset, 10));

		olist.add(createHistogram("log odds normal adjusted > "
				+ KarkinosProp.LognThres, "odds score", "count", dataset, 12,
				KarkinosProp.LognThres));

		olist.add(createHistogram("log odds tumor", "odds score", "count",
				dataset, 11));

		olist.add(createHistogram("log odds tumor adjusted >"
				+ +KarkinosProp.LogtThres, "odds score", "count", dataset, 13,
				KarkinosProp.LogtThres));

		return olist;

	}

	private static JFreeChart createHistogram(String title, String xlabel,
			String ylabel, DataSet dataset, int flg) {

		HistogramDataset hdataset = getDataSet(dataset, flg);
		JFreeChart chart = ChartFactory.createHistogram(title, xlabel, ylabel,
				hdataset, PlotOrientation.VERTICAL, true, false, false);
		Plot plot = chart.getPlot();
		plot.setBackgroundPaint(Color.WHITE);
		setSeriesPaint(chart);
		return chart;

	}

	private static void setSeriesPaint(JFreeChart chart) {
		try {
			XYItemRenderer xyr = ((XYPlot)chart.getPlot()).getRenderer();
			xyr.setSeriesPaint(0, Color.RED.darker());
			xyr.setSeriesPaint(1, Color.BLUE.darker());
			xyr.setSeriesPaint(2, Color.GREEN.darker());
			xyr.setSeriesPaint(3, Color.LIGHT_GRAY.darker());
			xyr.setSeriesPaint(4, Color.CYAN.darker());
		} catch (Exception ex) {
		}
		try {
			BarRenderer br = (BarRenderer)
							(((XYPlot)chart.getPlot()).getRenderer());
			br.setShadowVisible(false);
			br.setBarPainter(new StandardBarPainter());
			
		} catch (Exception ex) {
		}
		
	}

	private static JFreeChart createHistogram(String title, String xlabel,
			String ylabel, DataSet dataset, int flg, double line) {

		HistogramDataset hdataset = getDataSet(dataset, flg);
		JFreeChart chart = ChartFactory.createHistogram(title, xlabel, ylabel,
				hdataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		addMaker(plot, line);
		plot.setBackgroundPaint(Color.WHITE);
		setSeriesPaint(chart);
		
		try {
			plot.setRangeGridlinesVisible(true);
			plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
			plot.setDomainGridlinesVisible(true);
			plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
			BarRenderer br = (BarRenderer) plot.getRenderer();
			br.setShadowVisible(false);
			br.setBarPainter(new StandardBarPainter());

		} catch (Exception ex) {
		}
		return chart;

	}

	private static void addMaker(XYPlot xyplot, double d) {
		Marker marker0 = new ValueMarker(d);
		marker0.setPaint(Color.RED);
		xyplot.addDomainMarker(marker0);

	}

	private static double getData(SNVHolder snv, float tumorRratio, int flg) {
		//
		PileUPResult pir = snv.getTumor();
		boolean isindel = pir.isIndel();
		if (flg == 1) {

			float adjustedratio = CalcUtils.getTumorrateAdjustedRatio(
					snv, tumorRratio);
			return adjustedratio;

		} else if (flg == 2) {
			return pir.getRatio();
		} else if (flg == 3) {
			if (!isindel) {
				return pir.getBQrms();
			}

		} else if (flg == 4) {
			return pir.getMQrms();
		} else if (flg == 5) {
			return pir.getPhred();
		} else if (flg == 6) {
			return snv.getFilterResult().getSeqEntropy();
		} else if (flg == 7) {
			return snv.getFilterResult().getMappability();
		} else if (flg == 8) {
			double d = snv.getFilterResult().getPval4FiserDirectional();
			if (d == 0)
				d = -1;
			return d;
		} else if (flg == 9) {
			return snv.getPvalFisher();
		} else if (flg == 10) {
			if (!isindel) {
				return snv.getFilterResult().getLogn();
			}
		} else if (flg == 11) {
			if (!isindel) {
				return snv.getFilterResult().getLogt();
			}
		} else if (flg == 12) {
			if (!isindel) {
				return snv.getFilterResult().getLognAjusted();
			}
		} else if (flg == 13) {
			if (!isindel) {
				return snv.getFilterResult().getLogtAjusted();
			}
		}

		return -1000;
	}

	private static HistogramDataset getDataSet(DataSet dataset, int flg) {

		HistogramDataset histdata = new HistogramDataset();
		float tumorRratio = dataset.getTumorRatio();

		List<Double> all = new ArrayList<Double>();
		List<Double> dbSNP = new ArrayList<Double>();
		List<Double> filterOK = new ArrayList<Double>();
		List<Double> filter2OK = new ArrayList<Double>();

		double min = 0;
		double max = 0;
		for (SNVHolder snv : dataset.getSnvlist()) {

			if (snv.getFilterResult() != null) {
				double d = getData(snv, tumorRratio, flg);
				if (d == -1000) {
					continue;
				}
				if (max == 0)
					max = d;
				if (min == 0)
					min = d;
				if (min > d) {
					min = d;
				}
				if (max < d) {
					max = d;
				}
				all.add(d);
				if (snv.getFilterResult().getDbSNPbean() != null) {
					dbSNP.add(d);
				}
				if (snv.getFilterResult().isPassFilter()) {
					filterOK.add(d);
				}
				if (CalcUtils.pass2(snv.getFilterResult())) {
					filter2OK.add(d);
				}
			}

		}
		if (min > 0)
			min = 0;
		//
		histdata.addSeries("final cand", toAry(filter2OK), 100, min, max);
		histdata.addSeries("candidate", toAry(filterOK), 100, min, max);
		histdata.addSeries("dbSNP", toAry(dbSNP), 100, min, max);
		histdata.addSeries("all", toAry(all), 100, min, max);
		return histdata;

	}

	private static double[] toAry(List<Double> list) {

		double[] ary = new double[list.size()];
		int idx = 0;
		for (double d : list) {

			ary[idx] = d;
			idx++;
		}
		return ary;
	}

}
