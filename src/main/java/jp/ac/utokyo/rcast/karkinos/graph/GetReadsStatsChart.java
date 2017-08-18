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
import java.io.FileInputStream;
import java.io.IOException;
import java.text.DateFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CombinedDomainCategoryPlot;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.CombinedRangeCategoryPlot;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.lowagie.text.BadElementException;
import com.lowagie.text.Cell;
import com.lowagie.text.Table;

import jp.ac.utokyo.rcast.karkinos.cmd.KarkinosCmd;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.DepthCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.readssummary.XYSeriesForReadsDepth;

public class GetReadsStatsChart {

	private final static int SIZE = 2;

	public static List<DisplayObject> getChartLists(ReadsSummary readsSummary,
			String readsStat) {

		//
		List<DisplayObject> list = new ArrayList<DisplayObject>();
		// Table all
		List olist1 = new ArrayList();
		addFileInfoStr(olist1, readsSummary);

		olist1.add(getSummaryTable(readsSummary, readsStat));
		// if (readsSummary.isPairStats()) {

		// }
		olist1.add(getCaptureCoveageTable(readsSummary));
		olist1.add(getCaptureCoveageGraph(readsSummary));
		olist1.add(getInsertSizeGraph(readsSummary));
		// olist1.add(getLibrarySizeTable());

		list.add(new DisplayObject(olist1, SIZE, "Reads stats"));

		// Table chr normal
		list.add(new DisplayObject(getChrSummaryTable(
				readsSummary.getNormalCounter(),
				readsSummary.getNormalPerChrom()), SIZE, "Normal Reads stats"));

		// Graph normal
		Object objnormal = getSummaryGraph(readsSummary.getNormalPerChrom(),
				"Normal");
		DisplayObject normalGraph = new DisplayObject(objnormal, SIZE,
				"Normal reads stats");
		list.add(normalGraph);

		// Table chr tumor
		list.add(new DisplayObject(
				getChrSummaryTable(readsSummary.getTumorCounter(),
						readsSummary.getTumorPerChrom()), SIZE,
				"Tumor Reads stats"));

		// Graph tumor
		Object objTumor = getSummaryGraph(readsSummary.getTumorPerChrom(),
				"Tumor");
		DisplayObject dobj2 = new DisplayObject(objTumor, SIZE,
				"Tumor reads stats");
		list.add(dobj2);
		return list;

	}

	private static JFreeChart getInsertSizeGraph(ReadsSummary readsSummary) {

		// make a common vertical axis for all the sub-plots
		NumberAxis xAxis = new NumberAxis("Library size");
		int[] interval = readsSummary.getInsertSizeInterval();
		xAxis.setRange(interval[0], interval[1]);
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis("counts");
		yAxis.setAutoRangeIncludesZero(false);

		final XYLineAndShapeRenderer renderer0 = new XYLineAndShapeRenderer();
		renderer0.setSeriesShapesVisible(0, false);
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setSeriesShapesVisible(1, false);
		renderer0.setSeriesPaint(1, ChartColor.RED);

		if (readsSummary.getNormalCounter().isPairStats()) {
			XYSeriesCollection data0 = new XYSeriesCollection();
			data0.addSeries(getXYSeries(readsSummary.getNormalCounter(),
					interval, "Normal")); // add subplot 1...
			data0.addSeries(getXYSeries(readsSummary.getTumorCounter(),
					interval, "Tumor"));
			XYPlot subplot1 = new XYPlot(data0, xAxis, yAxis, renderer0);
			subplot1.setDomainCrosshairVisible(true);
			subplot1.setRangeCrosshairVisible(true);
			JFreeChart chart = new JFreeChart("Library size",
					JFreeChart.DEFAULT_TITLE_FONT, subplot1, true);
			return chart;

		}

		return null;

	}

	private static XYSeries getXYSeries(ReadsCounter readsCounter,
			int[] interval, String label) {

		XYSeries series = new XYSeries(label);
		TreeMap<Integer, CounterA> inserSizeMap = readsCounter
				.getInserSizeMap();
		for (int n = interval[0]; n < interval[1]; n++) {
			int x = n;
			int cnt = 0;
			if (inserSizeMap.containsKey(x)) {
				cnt = inserSizeMap.get(x).getCnt();
			}
			series.add(x, cnt);
		}
		return series;
	}

	private static Object getCaptureCoveageGraph(ReadsSummary readsSummary) {

		// CategoryAxis domainAxis = new CategoryAxis("depth");
		// NumberAxis yAxis = new NumberAxis("counts");
		// yAxis.setAutoRangeIncludesZero(true);
		XYSeriesForReadsDepth depthSummary = new XYSeriesForReadsDepth(
				readsSummary);
		//
		// // make a horizontally combined plot
		// // make a horizontally combined plot
		// CombinedDomainCategoryPlot parent = new
		// CombinedDomainCategoryPlot(domainAxis);

		final JFreeChart chart = ChartFactory.createBarChart(
				"reads depth coveage", // chart title
				"depth", // domain axis label
				"count", // range axis label
				depthSummary.getXYSeries(), // data
				PlotOrientation.VERTICAL, // orientation
				true, // include legend
				true, // tooltips?
				false // URLs?
				);
		setNo3D(chart);

		// final JFreeChart chartTumor = ChartFactory.createBarChart(
		// "reads depth coveage(tumor)", // chart title
		// "depth", // domain axis label
		// "count", // range axis label
		// depthSummary.getTumorXYSeries(), // data
		// PlotOrientation.VERTICAL, // orientation
		// true, // include legend
		// true, // tooltips?
		// false // URLs?
		// );
		CategoryPlot normalPlot = (CategoryPlot) chart.getPlot();
		normalPlot.setBackgroundPaint(Color.WHITE);
		// normalPlot.setDomainAxis(domainAxis);

		// CategoryPlot tumorPlot = (CategoryPlot)chartTumor.getPlot();
		// tumorPlot.setBackgroundPaint(Color.WHITE);
		// tumorPlot.setDomainAxis(domainAxis);

		// parent.add(normalPlot, 1);
		// parent.add(tumorPlot, 1);
		// createLineChart(data, string, string2, string3);
		// JFreeChart chart = new JFreeChart("InsertSize",
		// JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;

	}

	private static float getX20percent(DepthCounter dc) {

		try {
			long total = dc.getTotal();
			Map<Integer, CounterA> mt = dc.getMap();
			int[] keys = new int[] { 0, 1, 10 };
			long less20x = 0;
			for (int n : keys) {
				CounterA ca = mt.get(n);
				if (ca != null) {
					less20x = ca.getCnt();
				}

			}
			return (float) (((double) (total-less20x) / (double) total) * 100);
		} catch (Exception ex) {

		}
		return 0f;
	}

	private static Object getCaptureCoveageTable(ReadsSummary readsSummary) {
		Table table;
		try {
			table = new Table(5);

			// row1
			table.addCell("coveage");
			table.addCell("normal count");
			table.addCell("normal %");
			table.addCell("tumor count");
			table.addCell("tumor %");

			DepthCounter nd = readsSummary.getNormalDepth();
			long normaltotal = nd.getTotal();
			Map<Integer, CounterA> mn = nd.getMap();
			DepthCounter td = readsSummary.getTumorDepth();
			long tumortotal = td.getTotal();
			Map<Integer, CounterA> mt = td.getMap();

			Set<Integer> keySet = new TreeSet<Integer>();
			keySet.addAll(mn.keySet());
			keySet.addAll(mt.keySet());

			for (Integer key : keySet) {

				int countNormal = 0;
				float normalP = 0;
				int countTumor = 0;
				float tumorP = 0;

				String s;
				int n = key;
				if (n >= 1000) {
					s = "1000 -";
					// } else if (n >= 300) {
					// s = "" + n + "-" + (n + 99);
				} else if (n >= 10) {
					s = "" + n + "-" + (n + 9);
				} else if (n >= 1) {
					s = "" + n + "-" + (n + 8);
				} else {
					s = "" + n;
				}
				CounterA normal = mn.get(key);
				if (normal != null) {
					countNormal = normal.getCnt();
				}
				CounterA tumor = mt.get(key);
				if (tumor != null) {
					countTumor = tumor.getCnt();
				}
				normalP = (float) (((double) countNormal / (double) normaltotal) * 100);
				tumorP = (float) (((double) countTumor / (double) tumortotal) * 100);

				table.addCell(s);

				table.addCell(format(countNormal));
				table.addCell(format(normalP) + "%");

				table.addCell(format(countTumor));
				table.addCell(format(tumorP) + "%");

			}

			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	private static void addFileInfoStr(List olist1, ReadsSummary readsSummary) {

		olist1.add("");
		olist1.add(KarkinosCmd.getVersion().orElse(""));
		Date date = new Date();
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss.SSS");
		olist1.add(df.format(date));
		olist1.add("");
		olist1.add("normal: " + readsSummary.getNormalbam());
		olist1.add("tumor: " + readsSummary.getTumorbam());
		olist1.add("target: " + readsSummary.getTaretbed());
		olist1.add("");
		// olist1.add("normal mean library size: " +
		// readsSummary.getNormalCounter().getMeanInsertSize());
		// olist1.add("tumor mean library size: " +
		// readsSummary.getTumorCounter().getMeanInsertSize());

		olist1.add("used properties");
		olist1.add("");
		olist1.add(KarkinosProp.getInfoString());

	}

	private static Object getSummaryGraph(Map<String, ReadsCounter> datamap,
			String title) {

		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (Entry<String, ReadsCounter> entry : datamap.entrySet()) {

			String chrom = entry.getKey();
			ReadsCounter rc = entry.getValue();
			int target = rc.getTotalOnTarget();
			int total = rc.getTotalmap();
			int outoftarget = total - target;
			dataset.addValue(target, "on target", chrom);
			dataset.addValue(outoftarget, "out of target", chrom);

		}

		JFreeChart jfc = ChartFactory.createStackedBarChart(title, "chrom",
				"tag counts", dataset, PlotOrientation.VERTICAL, true, false,
				false);

		jfc.getPlot().setBackgroundPaint(Color.WHITE);
		setNo3D(jfc);
		return jfc;
	}

	private static Object getLibrarySizeTable(ReadsSummary readsSummary) {

		Table table;
		try {
			table = new Table(6);

			// row1
			table.addCell("Library size");
			table.addCell("normal");
			table.addCell("tumor");

			table.addCell("");
			table.addCell(format(readsSummary.getNormalCounter()
					.getMeanInsertSize()));
			table.addCell(format(readsSummary.getTumorCounter()
					.getMeanInsertSize()));
			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	private static Object getChrSummaryTable(ReadsCounter totalc,
			Map<String, ReadsCounter> datamap) {

		Table table;
		try {
			table = new Table(6);

			// row1
			table.addCell("chr/");
			table.addCell("on target");
			table.addCell("out of target");
			table.addCell("on tag unique");
			table.addCell("outoftag unique");
			table.addCell("ontag %");

			for (Entry<String, ReadsCounter> entry : datamap.entrySet()) {

				String chrom = entry.getKey();
				//
				boolean usualchrom = usualChrom(chrom);
				
				ReadsCounter rc = entry.getValue();
				int target = rc.getTotalOnTarget();
				int total = rc.getTotalmap();
				int outoftarget = total - target;
				int totalu = rc.getTotalunique();
				int uniqueoutoftarget = totalu - rc.getTotalUniqueOntarget();

				table.addCell(chrom);
				table.addCell(format(target));
				table.addCell(format(outoftarget));
				table.addCell(format(rc.getTotalUniqueOntarget()));
				table.addCell(format(uniqueoutoftarget));
				table.addCell(format(rc.onTargetParcent()) + "%");

			}

			// write total
			table.addCell("Total");
			ReadsCounter rc = totalc;
			int target = rc.getTotalOnTarget();
			int total = rc.getTotalmap();
			int outoftarget = total - target;
			int totalu = rc.getTotalunique();
			int uniqueoutoftarget = totalu - rc.getTotalUniqueOntarget();
			table.addCell(format(target));
			table.addCell(format(outoftarget));
			table.addCell(format(rc.getTotalUniqueOntarget()));
			table.addCell(format(uniqueoutoftarget));
			table.addCell(format(rc.onTargetParcent()) + "%");

			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}
	
	private static boolean usualChrom(String chrom) {

		String chromnum = chrom;
		if (chromnum.contains("chr")) {
			chromnum = chromnum.replaceAll("chr", "");
		}
		if (StringUtils.isNumeric(chromnum))
			return true;
		if (chromnum.equalsIgnoreCase("X"))
			return true;
		if (chromnum.equalsIgnoreCase("Y"))
			return true;
		if (chromnum.equalsIgnoreCase("M"))
			return true;
		if (chromnum.equalsIgnoreCase("mt"))
			return true;

		return false;

	}

	private static Object getSummaryTable(ReadsSummary readsSummary,
			String readsStat) {

		try {

			ReadsCounter normal = readsSummary.getNormalCounter();
			ReadsCounter tumor = readsSummary.getTumorCounter();

			Table table = new Table(3);
			int[] width = new int[] { 60, 20, 20 };
			table.setWidths(width);

			// row1
			table.addCell("reads stats/");
			table.addCell("normal");
			table.addCell("tumor");

			// // row
			if (readsStat != null) {

				String[] sa = readsStat.split(",");
				String normalfile = sa[0];
				Properties np = getProp(normalfile);
				String tumorfile = sa[1];
				Properties tp = getProp(tumorfile);
				// datamap.put("total", total);
				// datamap.put("passfilter", passfilter);
				// datamap.put("mapped", mapped);
				// datamap.put("nonduplicate", nonduplicate);
				// datamap.put("duplicate", duplicate);
				// //datamap.put("identityLow", identityLow);
				// datamap.put("unique", unique);
				// datamap.put("proper", proper);
				// datamap.put("properOrunique", properOrunique);
				// datamap.put("bincount", 1l);
				// datamap.put("afterrealgin", afterrealgin);
				// datamap.put("identityLow", identityLow);
				String[] keys = new String[] { "total", "passfilter", "mapped",
						"nonduplicate", "duplicate", "unique", "proper",
						"properOrunique", "afterrealgin", "identityLow" };
				for (String key : keys) {
					if (np.containsKey(key)) {

						// row2
						String kyeDisp = key.replaceAll("Or", " or ");
						table.addCell(kyeDisp + " tags");
						table.addCell(format(np.getProperty(key)));
						table.addCell(format(tp.getProperty(key)));

						if (key.equals("duplicate")) {

							try {
								// row 13
								table.addCell("duplicate (%)");
								table.addCell(format(getParcent(
										np.getProperty(key),
										np.getProperty("nonduplicate")))
										+ "%");
								table.addCell(format(getParcent(
										tp.getProperty(key),
										tp.getProperty("nonduplicate")))
										+ "%");
							} catch (Exception ex) {
							}
						}
					}
				}

			}

			// row
			table.addCell("Total used tags");
			table.addCell(format(normal.getTotalmap()));
			table.addCell(format(tumor.getTotalmap()));

			// row4
			table.addCell("Target tags");
			table.addCell(format(normal.getTotalOnTarget()));
			table.addCell(format(tumor.getTotalOnTarget()));
			// row5
			table.addCell("Target unique tags");
			table.addCell(format(normal.getTotalUniqueOntarget()));
			table.addCell(format(tumor.getTotalUniqueOntarget()));

			if (readsSummary.isPairStats()) {

				// row6
				table.addCell("First reads");
				table.addCell(format(normal.getFirstReads()));
				table.addCell(format(tumor.getFirstReads()));
				// row7
				table.addCell("Second reads");
				table.addCell(format(normal.getSecondReads()));
				table.addCell(format(tumor.getSecondReads()));
				// row8
				table.addCell("Both mapped");
				table.addCell(format(normal.getBothmap()));
				table.addCell(format(tumor.getBothmap()));
				// row9
				table.addCell("Proper reads");
				table.addCell(format(normal.getPropermap()));
				table.addCell(format(tumor.getPropermap()));
				// row10
				table.addCell("mean library size");
				table.addCell(format(normal.getMeanInsertSize()));
				table.addCell(format(tumor.getMeanInsertSize()));

			}

			// // row11
			// row 12
			table.addCell("mean depth");
			table.addCell(format(readsSummary.getNormalDepth().getMeanDepth()));
			table.addCell(format(readsSummary.getTumorDepth().getMeanDepth()));
			
			table.addCell("total covered region");
			table.addCell(format(readsSummary.getNormalDepth().getTotal()));
			table.addCell(format(readsSummary.getTumorDepth().getTotal()));
			
			// row 13
			table.addCell("On target(%)");
			table.addCell(format(normal.onTargetParcent()) + "%");
			table.addCell(format(tumor.onTargetParcent()) + "%");

			table.addCell("more than X20(%)");
			table.addCell(format(getX20percent(readsSummary.getNormalDepth()))
					+ "%");
			table.addCell(format(getX20percent(readsSummary.getTumorDepth()))
					+ "%");
			
			
			table.addCell("mean depth CDS");
			table.addCell(format(readsSummary.getNormalDepth().getMeanCDSDepth()));
			table.addCell(format(readsSummary.getTumorDepth().getMeanCDSDepth()));
			
			table.addCell("more than X10(%) CDS");
			table.addCell(format(readsSummary.getNormalDepth().getOver10X())+ "%");
			table.addCell(format(readsSummary.getTumorDepth().getOver10X())+ "%");
			
			table.addCell("more than X20(%) CDS");
			table.addCell(format(readsSummary.getNormalDepth().getOver20X())+ "%");
			table.addCell(format(readsSummary.getTumorDepth().getOver20X())+ "%");
			

			
			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	private static double getParcent(String dup, String nondup) {

		double d = Double.parseDouble(dup);
		double nond = Double.parseDouble(nondup);
		d = (d / (d + nond)) * 100;
		return d;
	}

	private static Properties getProp(String file) throws IOException {
		// TODO Auto-generated method stub
		FileInputStream fis = new FileInputStream(file);
		Properties p = new Properties();
		p.load(fis);
		return p;
	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static String format(String num) {
		try {
			long l = Long.parseLong(num);
			return format(l);
		} catch (Exception ex) {

		}
		try {
			double d = Double.parseDouble(num);
			return format(d);
		} catch (Exception ex) {

		}
		return "";
	}

	private static String format(long num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static String format(int num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static void setNo3D(JFreeChart chart) {

		try {
			CategoryPlot cp = chart.getCategoryPlot();
			BarRenderer br = (BarRenderer) cp.getRenderer();
			br.setShadowVisible(false);
			br.setBarPainter(new StandardBarPainter());

		} catch (Exception ex) {
		}
		try {
			CategoryPlot cp = chart.getCategoryPlot();
			CategoryAxisSkipLabels ca = new CategoryAxisSkipLabels();
			cp.setDomainAxis(ca);

		} catch (Exception ex) {
		}
		//

	}
}
