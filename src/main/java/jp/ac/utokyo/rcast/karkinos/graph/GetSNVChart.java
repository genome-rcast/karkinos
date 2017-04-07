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

import java.awt.BasicStroke;
import java.awt.Color;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.lowagie.text.BadElementException;
import com.lowagie.text.Cell;
import com.lowagie.text.Paragraph;
import com.lowagie.text.Table;

import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.distribution.DataHolderByCN;
import jp.ac.utokyo.rcast.karkinos.distribution.XYDepthRatioExtract;
import jp.ac.utokyo.rcast.karkinos.distribution.XYSeriesExtract;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.graph.output.FormatHelper;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.MutationCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.readssummary.SNPDepthCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.SummaryStatsHolder;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.GenotypeKeyUtils;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.*;

public class GetSNVChart {

	private final static int SIZE = 1; 

	public static List<DisplayObject> getChartLists(DataSet dataset, PeaksInfo pi) {

		AnalyseDist analyseDist = dataset.getAnalyseDist();
		Map<Float, DataHolderByCN> snvMAP = analyseDist.getMap();
		XYSeriesExtract xyse = new XYSeriesExtract(snvMAP);
		XYDepthRatioExtract xyrd = new XYDepthRatioExtract(dataset);
		
		List<DisplayObject> list = new ArrayList<DisplayObject>();

		list.add(new DisplayObject(getTRTable(dataset,pi), SIZE,
				"Estimated Tumor Contents"));

		List[] stats = getFilterStatsTables(dataset);
		
		list.add(new DisplayObject(stats[0], 2,
				"Filtered Candidates"));
		
		
		//
		List oList = new ArrayList();
		oList.add(getRatioDepthGraph(xyrd.getTumorRatioDepth(0)));
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.25f) {
			IntervalXYDataset normal = xyse.getXYSeries(degree, true, true,
					false);
			IntervalXYDataset tumor = xyse.getXYSeries(degree, true, false,
					false);
			
			if(normal==null||tumor==null){
				continue;
			}
			
			if (degree == 0.5f) {

				float observeRatio = dataset.getAnalyseDist()
						.getTumorratioFromLOH().getObservedratio();
				oList.add(getHetroSNPGraph(normal, tumor, degree, observeRatio));
			} else if(degree == 1.5f){
				float observeRatio = dataset.getAnalyseDist()
		 		.getTumorratioFromGAIN().getObservedratio();
				oList.add(getHetroSNPGraph(normal, tumor, degree, observeRatio));
			}else{	
				oList.add(getHetroSNPGraph(normal, tumor, degree));
			}
			
		}
		list.add(new DisplayObject(oList, SIZE, "Hetro SNP Graph"));

		List oList2 = new ArrayList();
		oList2.add(getRatioDepthGraph(xyrd.getTumorRatioDepth(1)));
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.25f) {
			IntervalXYDataset normal = xyse.getXYSeries(degree, false, true,
					false);
			IntervalXYDataset tumor = xyse.getXYSeries(degree, false, false,
					false);

			if(normal==null||tumor==null){
				continue;
			}

			oList2.add(getHetroSNPGraph(normal, tumor, degree));
		}
		list.add(new DisplayObject(oList2, SIZE, "mutation Graph"));

		List oList3 = new ArrayList();
		oList3.add(getRatioDepthGraph(xyrd.getTumorRatioDepth(2)));
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.25f) {
			IntervalXYDataset normal = xyse.getXYSeries(degree, false, true,
					true);
			IntervalXYDataset tumor = xyse.getXYSeries(degree, false, false,
					true);
			if(normal==null||tumor==null){
				continue;
			}

			oList3.add(getHetroSNPGraph(normal, tumor, degree));
		}
		list.add(new DisplayObject(oList3, SIZE,
				"mutation Graph after filtering"));

		List oList4 = new ArrayList();
		oList4.add(getRatioDepthGraph(xyrd.getTumorRatioDepth(3)));
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.25f) {
			IntervalXYDataset normal = xyse.getXYSeriesAfterFinalFilter(degree,
					true);
			IntervalXYDataset tumor = xyse.getXYSeriesAfterFinalFilter(degree,
					false);
			if(normal==null||tumor==null){
				continue;
			}			
			oList4.add(getHetroSNPGraph(normal, tumor, degree));
		}
		list.add(new DisplayObject(oList4, SIZE,
				"mutation Graph after bayesian filtering"));

		list.add(new DisplayObject(stats[1], 1,"genotype info"));

		
		return list;
	}

	private static JFreeChart getRatioDepthGraph(XYSeriesCollection data) {
		
		JFreeChart chart 
			= ChartFactory.createScatterPlot("ratio/depth",
					"ratio", "read depth" ,data, 
					PlotOrientation.VERTICAL, true, false, false);
		
		XYPlot plot = chart.getXYPlot();
		plot.getDomainAxis().setRange(0, 1);
		plot.setBackgroundPaint(Color.WHITE);
		plot.setRangeGridlinesVisible(true);
		plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
		plot.setDomainGridlinesVisible(true);
		plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
		return chart;
		
	}

	static final int PASS_FILTERFinal = 1000;


	

	private static List[] getFilterStatsTables(DataSet dataset) {

		List olist = new ArrayList();
		List olistG = new ArrayList();
		Map<Integer, CounterA> sm = new LinkedHashMap<Integer, CounterA>();
		Map<Integer, CounterA> ti = new LinkedHashMap<Integer, CounterA>();

		Map<Integer, CounterA> smsolo = new LinkedHashMap<Integer, CounterA>();
		Map<Integer, CounterA> tisolo = new LinkedHashMap<Integer, CounterA>();

		//
		int insersion = 0;
		int delation = 0;
		int totalmutaion = 0;
		//
		Map<String, MutationCounter> mutationCounter = new TreeMap<String, MutationCounter>();
		Map<String, SNPDepthCounter> mutationCounterNormalSNP = new TreeMap<String, SNPDepthCounter>();
		Map<String, SNPDepthCounter> normalmutationCounter = new TreeMap<String, SNPDepthCounter>();
		Map<String,CounterA> chrCounts = new LinkedHashMap<String,CounterA>();
		
		//		String keys[] = new String[] { "AtoT", "AtoG", "AtoC", "TtoA", "TtoG",
//				"TtoC", "CtoA", "CtoT", "CtoG", "GtoA", "GtoT", "GtoC" };

		
		for (String key : GenotypeKeyUtils.keys1) {
			
			mutationCounter.put(key, new MutationCounter());
			mutationCounterNormalSNP.put(key, new SNPDepthCounter());
			normalmutationCounter.put(key, new SNPDepthCounter());
		}

		int mutationfilterout = 0;
		int indelfilterout = 0;

		for (SNVHolder snv : dataset.getSnvlist()) {

			int flg = snv.getFlg();
			if (flg == PileUP.NormalSNP || flg == PileUP.REGBOTH) {
				char ref = snv.getNormal().getGenomeR();
				char alt = snv.getNormal().getALT();
				String key = ref + "to" + alt;
				key = GenotypeKeyUtils.aggrigateKeys(key);
				SNPDepthCounter counter = null;
				if (normalmutationCounter.containsKey(key)) {
					counter = normalmutationCounter.get(key);
					counter.reg(snv.getNormal().getTotalcnt());
				}
			}
			if (flg == PileUP.NormalSNP) {

				//
				
				char ref = snv.getNormal().getGenomeR();
				char alt = snv.getNormal().getALT();
				String key = ref + "to" + alt;
				key = GenotypeKeyUtils.aggrigateKeys(key);
				SNPDepthCounter counter = null;
				boolean hetroSNP = snv.getNormal().getRatio() < 0.6;
				if (hetroSNP && mutationCounterNormalSNP.containsKey(key)) {
					counter = mutationCounterNormalSNP.get(key);
					counter.reg(snv.getNormal().getTotalcnt());
				}

			} else if ((flg == PileUP.SomaticMutation)
					|| (flg == PileUP.TumorINDEL)) {

				if (null != snv.getFilterResult()) {
					
					
					if(snv.getChr().equals("chr1")){
						if(snv.getPos()==13645209){
							System.out.println("here");
						}
					}
					Set<Integer> filterflgs = snv.getFilterResult()
							.getPassFilterFlg();
					Set<Integer> infoflgs = snv.getFilterResult().getInfoFlg();
						
					Map<Integer, CounterA> map = (flg == PileUP.SomaticMutation) ? sm
							: ti;

					Map<Integer, CounterA> mapsolo = (flg == PileUP.SomaticMutation) ? smsolo
							: tisolo;

					CounterA ca = null;
					if (snv.getFilterResult().isPassFilter()) {

						if (map.containsKey(PASS_FILTER)) {
							ca = map.get(PASS_FILTER);
						} else {
							ca = new CounterA();
							map.put(PASS_FILTER, ca);
						}
						ca.inc();

						if (CalcUtils.pass2(snv.getFilterResult())) {
							
							addchrCounts(chrCounts,snv.getChr());
							//
							CounterA caf = null;
							if (map.containsKey(PASS_FILTERFinal)) {
								caf = map.get(PASS_FILTERFinal);
							} else {
								caf = new CounterA();
								map.put(PASS_FILTERFinal, caf);
							}
							caf.inc();
						}else{
							
							for (int filterflg : infoflgs) {
								if (map.containsKey(filterflg)) {
									ca = map.get(filterflg);
								} else {
									ca = new CounterA();
									map.put(filterflg, ca);
								}
								ca.inc();
							}

							// filtered by this property only
							if (infoflgs.size() == 1) {
								int filterflg = infoflgs.iterator().next();
								CounterA counter2 = null;
								if (mapsolo.containsKey(filterflg)) {
									counter2 = mapsolo.get(filterflg);
								} else {
									counter2 = new CounterA();
									mapsolo.put(filterflg, counter2);
								}
								counter2.inc();
							}							
							
						}

						if (flg == PileUP.TumorINDEL) {
							if (snv.getTumor().isInsersion()) {
								insersion++;
							} else {
								delation++;
							}

						} else {
							char ref = snv.getTumor().getGenomeR();
							char alt = snv.getTumor().getALT();
							String key = ref + "to" + alt;
							key = GenotypeKeyUtils.aggrigateKeys(key);
							MutationCounter counter = null;
							if (mutationCounter.containsKey(key)) {
								counter = mutationCounter.get(key);
								int copynumber = (int) (snv.getCi()
										.getVaridateVal() * 2);
								PileUPResult pir = snv.getTumor();
								float ratio = pir.getRatio();
								float fpval = (float) snv.getPvalFisher();
								float normallogodd = snv.getNormal()
										.getRefLogLikeHood();
								float adgustedmutateodds = (float) snv
										.getFilterResult().getLogtAjusted();
								if (fpval == 0) {
									fpval = 0.0000000000000001f;
								}
								float logp = (float) (-1 * Math.log10(fpval));
								int depth = pir.getTotalcnt();
								counter.regStat(copynumber, ratio, logp, depth);
								if (CalcUtils.pass2(snv.getFilterResult())) {
									counter.regStat2(depth);
								}
								totalmutaion++;
							}
							// /

						}

					} else {
						// filter out
						if (flg == PileUP.TumorINDEL) {
							indelfilterout++;
						} else {
							mutationfilterout++;
						}
						for (int filterflg : filterflgs) {
							if (map.containsKey(filterflg)) {
								ca = map.get(filterflg);
							} else {
								ca = new CounterA();
								map.put(filterflg, ca);
							}
							ca.inc();
						}

						// filtered by this property only
						if (filterflgs.size() == 1) {
							int filterflg = filterflgs.iterator().next();
							CounterA counter2 = null;
							if (mapsolo.containsKey(filterflg)) {
								counter2 = mapsolo.get(filterflg);
							} else {
								counter2 = new CounterA();
								mapsolo.put(filterflg, counter2);
							}
							counter2.inc();
						}

					}
				}

			}
		}

		try {

			Table table1 = summarytable(sm, ti, smsolo, tisolo,
					mutationfilterout, indelfilterout);
			
			olist.add(table1);
			//
			
			// //////////////////
			// mutation stat
			Table tablemutation = new Table(3);
			tablemutation.addCell("mutation type");
			tablemutation.addCell("counts");
			tablemutation.addCell("%");

			Set<Entry<String, MutationCounter>> s = mutationCounter.entrySet();
			for (Entry<String, MutationCounter> e : s) {

				tablemutation.addCell(GenotypeKeyUtils.toDispKey(e.getKey()));
				MutationCounter mc = e.getValue();
				int cnt = mc.getTotal();
				float parcent = (float) ((double) cnt * 100 / (double) totalmutaion);
				tablemutation.addCell(String.valueOf(cnt));
				tablemutation.addCell(String.valueOf(parcent) + "%");

			}
	
			
			
			tablemutation.addCell("delation");
			tablemutation.addCell(String.valueOf(delation));
			tablemutation.addCell("-");

			tablemutation.addCell("insertion");
			tablemutation.addCell(String.valueOf(insersion));
			tablemutation.addCell("-");

			olist.add(tablemutation);
			
			
			///chr table
			//
			Table tablechr = new Table(2);
			tablechr.addCell("chr");
			tablechr.addCell("final candidates");
			Set<Entry<String, CounterA>> chrs = chrCounts.entrySet();
			for (Entry<String, CounterA> e : chrs) {

				tablechr.addCell(e.getKey());				
				tablechr.addCell(""+e.getValue().getCnt());

			}		
			
			olist.add(tablechr);

			// mutation stat
			Table tableallele = new Table(13);
			tableallele.addCell("type");
			tableallele.addCell("counts(n=1)");
			tableallele.addCell("ratio mean(n=1)");
			tableallele.addCell("ratio sd(n=1)");
			tableallele.addCell("counts(n=2)");
			tableallele.addCell("ratio mean(n=2)");
			tableallele.addCell("ratio sd(n=2)");
			tableallele.addCell("counts(n=3)");
			tableallele.addCell("ratio mean(n=3)");
			tableallele.addCell("ratio sd(n=3)");
			tableallele.addCell("counts(n=4)");
			tableallele.addCell("ratio mean(n=4)");
			tableallele.addCell("ratio sd(n=4)");
			//
			for (Entry<String, MutationCounter> e : s) {

				MutationCounter mc = e.getValue();
				// int cnt = mc.getTotal();
				tableallele.addCell(GenotypeKeyUtils.toDispKey(e.getKey()));
				for (int n = 1; n <= 4; n++) {
					SummaryStatsHolder ssh = mc.getSummaryStatsHolder(n);
					if (null != ssh) {
						tableallele.addCell(ssh.getRatio().getN() + "");
						tableallele.addCell(format(ssh.getRatio().getMean()));
						tableallele.addCell(format(ssh.getRatio()
								.getStandardDeviation()));
					} else {
						tableallele.addCell("-");
						tableallele.addCell("-");
						tableallele.addCell("-");
					}
				}

			}

			olist.add(tableallele);

			// mutation stat
			Table tableallelep = new Table(13);
			tableallelep.addCell("type");
			tableallelep.addCell("counts(n=1)");
			tableallelep.addCell("log(p val) mean(n=1)");
			tableallelep.addCell("log(p val) sd(n=1)");
			tableallelep.addCell("counts(n=2)");
			tableallelep.addCell("log(p val) mean(n=2)");
			tableallelep.addCell("log(p val) sd(n=2)");
			tableallelep.addCell("counts(n=3)");
			tableallelep.addCell("log(p val) mean(n=3)");
			tableallelep.addCell("log(p val) sd(n=3)");
			tableallelep.addCell("counts(n=4)");
			tableallelep.addCell("log(p val) mean(n=4)");
			tableallelep.addCell("log(p val) sd(n=4)");
			//
			for (Entry<String, MutationCounter> e : s) {

				MutationCounter mc = e.getValue();
				// int cnt = mc.getTotal();
				tableallelep.addCell(GenotypeKeyUtils.toDispKey(e.getKey()));
				for (int n = 1; n <= 4; n++) {
					SummaryStatsHolder ssh = mc.getSummaryStatsHolder(n);
					if (null != ssh) {

						tableallelep.addCell(ssh.getOddsratio().getN() + "");
						tableallelep.addCell(format(ssh.getOddsratio()
								.getMean()));
						tableallelep.addCell(format(ssh.getOddsratio()
								.getStandardDeviation()));

					} else {
						tableallelep.addCell("-");
						tableallelep.addCell("-");
						tableallelep.addCell("-");
					}
				}

			}

			olist.add(tableallelep);
			Set<Entry<String, SNPDepthCounter>> sNormal = normalmutationCounter 
					.entrySet();
			DefaultCategoryDataset cdatasetNormal = new DefaultCategoryDataset();

			for (Entry<String, SNPDepthCounter> e : sNormal) {

				SNPDepthCounter sc = e.getValue();
				int lowcnt = sc.getLowdepthcnt();
				cdatasetNormal.addValue(lowcnt, "x1-x2(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x2count = sc.getX2count();
				cdatasetNormal.addValue(x2count, "x2-x3(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x3count = sc.getX3count();
				cdatasetNormal.addValue(x3count, "x3-x4(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x4count = sc.getX2count();
				cdatasetNormal.addValue(x4count, "x4-x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int highdepth = sc.getHighdepthcnt();
				cdatasetNormal.addValue(highdepth, "< x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));

			}

			JFreeChart jfcNormal = ChartFactory.createStackedBarChart(
					"normal mutation frequency", "mutation", "counts", cdatasetNormal,
					PlotOrientation.VERTICAL, true, false, false);
			jfcNormal.getPlot().setBackgroundPaint(Color.WHITE);
			setSeriespaint(jfcNormal);
			
			olistG.add(jfcNormal);

			//
			Set<Entry<String, SNPDepthCounter>> sSNP = mutationCounterNormalSNP
					.entrySet();
			DefaultCategoryDataset cdatasetSNP = new DefaultCategoryDataset();

			for (Entry<String, SNPDepthCounter> e : sSNP) {

				SNPDepthCounter sc = e.getValue();
				int lowcnt = sc.getLowdepthcnt();
				cdatasetSNP.addValue(lowcnt, "x1-x2(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x2count = sc.getX2count();
				cdatasetSNP.addValue(x2count, "x2-x3(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x3count = sc.getX3count();
				cdatasetSNP.addValue(x3count, "x3-x4(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x4count = sc.getX2count();
				cdatasetSNP.addValue(x4count, "x4-x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int highdepth = sc.getHighdepthcnt();
				cdatasetSNP.addValue(highdepth, "< x5(" + KarkinosProp.mindepth
						+ ")",GenotypeKeyUtils.toDispKey(e.getKey()));

			}

			JFreeChart jfcSNP = ChartFactory.createStackedBarChart(
					"normal SNP frequency", "mutation", "counts", cdatasetSNP,
					PlotOrientation.VERTICAL, true, false, false);
			jfcSNP.getPlot().setBackgroundPaint(Color.WHITE);
			setSeriespaint(jfcSNP);

			olistG.add(jfcSNP);

			//
			DefaultCategoryDataset cdataset = new DefaultCategoryDataset();
			for (Entry<String, MutationCounter> e : s) {

				MutationCounter mc = e.getValue();
				int lowcnt = mc.getLowdepthcnt();
				cdataset.addValue(lowcnt, "x1-x2(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x2count = mc.getX2count();
				cdataset.addValue(x2count, "x2-x3(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x3count = mc.getX3count();
				cdataset.addValue(x3count, "x3-x4(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x4count = mc.getX2count();
				cdataset.addValue(x4count, "x4-x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int highdepth = mc.getHighdepthcnt();
				cdataset.addValue(highdepth, "< x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));

			}

			JFreeChart jfc = ChartFactory.createStackedBarChart(
					"mutation frequency", "mutation", "counts", cdataset,
					PlotOrientation.VERTICAL, true, false, false);
			jfc.getPlot().setBackgroundPaint(Color.WHITE);
			setSeriespaint(jfc);
			olistG.add(jfc);

			//
			DefaultCategoryDataset cdataset2 = new DefaultCategoryDataset();
			for (Entry<String, MutationCounter> e : s) {

				MutationCounter mc = e.getValue();
				int lowcnt = mc.getLowdepthcntH();
				cdataset2.addValue(lowcnt, "x1-x2(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x2count = mc.getX2countH();
				cdataset2.addValue(x2count, "x2-x3(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x3count = mc.getX3countH();
				cdataset2.addValue(x3count, "x3-x4(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int x4count = mc.getX2countH();
				cdataset2.addValue(x4count, "x4-x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));
				int highdepth = mc.getHighdepthcntH();
				cdataset2.addValue(highdepth, "< x5(" + KarkinosProp.mindepth
						+ ")", GenotypeKeyUtils.toDispKey(e.getKey()));

			}

			JFreeChart jfc2 = ChartFactory.createStackedBarChart(
					"mutation frequency filter2", "mutation", "counts",
					cdataset2, PlotOrientation.VERTICAL, true, false, false);
			jfc2.getPlot().setBackgroundPaint(Color.WHITE);
			setSeriespaint(jfc2);
			olistG.add(jfc2);

			return new List[]{olist,olistG};

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	private static Table summarytable(Map<Integer, CounterA> sm,
			Map<Integer, CounterA> ti, Map<Integer, CounterA> smsolo,
			Map<Integer, CounterA> tisolo, int mutationfilterout,
			int indelfilterout) throws BadElementException {
		Table table = new Table(3);
		//
		// row1
		table.addCell("filter type");
		table.addCell("somatic mutation");
		table.addCell("somatic indel");

		// public static final int SNP = 1;
		// public static final int Low_mappability = 5;
		// public static final int Low_complexty = 6;
		// public static final int Low_SNPQual = 7;
		// public static final int Basian_Filteling = 9;
		// public static final int SUPPORTED_READSNG = 19;
		// public static final int SUPPORTED_BY_ONEDirection = 11;
		// public static final int CONTAIN_Reccurent_MISMATCH = 12;
		// public static final int READSENDSONLY = 13;
		// public static final int TooManyMismatchReads = 14;
		// public static final int MutationAtSameCycle = 15;
		//
		table.addCell("Pass filter Filter2");
		table.addCell(getContStr(sm, PASS_FILTERFinal));
		table.addCell(getContStr(ti, PASS_FILTERFinal));

		table.addCell("Pass filter Filter1");
		table.addCell(getContStr(sm, FilterResult.PASS_FILTER));
		table.addCell(getContStr(ti, FilterResult.PASS_FILTER));

		table.addCell("Filtered out total");
		table.addCell(String.valueOf(mutationfilterout));
		table.addCell(String.valueOf(indelfilterout));
		
//		table.addCell("High normal TC adjuated ratio");
//		table.addCell(getContStr(sm, smsolo, FilterResult.High_normal_adjustedRatio));
//		table.addCell(getContStr(ti, tisolo, FilterResult.High_normal_adjustedRatio));
//		
		
		table.addCell("Bayesian filter for normal");
		table.addCell(getContStr(sm, smsolo, FilterResult.INFO_LOW_refOddsRatio));
		table.addCell(getContStr(ti, tisolo, FilterResult.INFO_LOW_refOddsRatio));
		
		table.addCell("Bayesian filter for tumor");
		table.addCell(getContStr(sm, smsolo, FilterResult.INFO_LOW_tumorOddsRatio));
		table.addCell(getContStr(ti, tisolo, FilterResult.INFO_LOW_tumorOddsRatio));
		
		table.addCell("Low adjuated allele ratio");
		table.addCell(getContStr(sm, smsolo, FilterResult.INFO_adjustAlleleFreq));
		table.addCell(getContStr(ti, tisolo, FilterResult.INFO_adjustAlleleFreq));
		
		table.addCell("Low support reads for TP rescued candidate");
		table.addCell(getContStr(sm, smsolo, FilterResult.INFO_adjustLowdepth));
		table.addCell(getContStr(ti, tisolo, FilterResult.INFO_adjustLowdepth));
		
		table.addCell("min support reads for non cosmic sites");
		table.addCell(getContStr(sm, smsolo, FilterResult.INFO_minimumSupportReads));
		table.addCell(getContStr(ti, tisolo, FilterResult.INFO_minimumSupportReads));
		
		table.addCell("Low TC adjuated ratio");
		table.addCell(getContStr(sm, smsolo, FilterResult.Low_tumor_adjustedRatio));
		table.addCell(getContStr(ti, tisolo, FilterResult.Low_tumor_adjustedRatio));
		
		table.addCell("Low AF in AF depth matrix");
		table.addCell(getContStr(sm, smsolo, FilterResult.Low_tumor_adjustedReads));
		table.addCell(getContStr(ti, tisolo, FilterResult.Low_tumor_adjustedReads));

		table.addCell("High low base quality reads ratio");
		table.addCell(getContStr(sm, smsolo, FilterResult.HighLowQualReads));
		table.addCell(getContStr(ti, tisolo, FilterResult.HighLowQualReads));
		
		table.addCell("dbSNP");
		table.addCell(getContStr(sm, smsolo, FilterResult.SNP));
		table.addCell(getContStr(ti, tisolo, FilterResult.SNP));

		//
		table.addCell("Supported by only one directional reads");
		table.addCell(getContStr(sm, smsolo, SUPPORTED_BY_ONEDirection));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.SUPPORTED_BY_ONEDirection));

		//
		table.addCell("Support reads contain another recurrent mismatch");
		table.addCell(getContStr(sm, smsolo,
				FilterResult.CONTAIN_Reccurent_MISMATCH));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.CONTAIN_Reccurent_MISMATCH));

		table.addCell("Support reads allelic inblance");
		table.addCell(getContStr(sm, smsolo,
				FilterResult.noStrandSpecific));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.noStrandSpecific));
		
		table.addCell("Support reads contain toomany mismatch");
		table.addCell(getContStr(sm, smsolo,
				FilterResult.TooManyMismatchReads));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.TooManyMismatchReads));

		//
		table.addCell("Low Mappability ");
		table.addCell(getContStr(sm, smsolo, FilterResult.Low_mappability));
		table.addCell(getContStr(ti, tisolo, FilterResult.Low_mappability));

		table.addCell("Low Map Quality ");
		table.addCell(getContStr(sm, smsolo, FilterResult.Low_MapQuality));
		table.addCell(getContStr(ti, tisolo, FilterResult.Low_MapQuality));

		table.addCell("Low sequence complexty");
		table.addCell(getContStr(sm, smsolo, FilterResult.Low_complexty));
		table.addCell(getContStr(ti, tisolo, FilterResult.Low_complexty));
		
		table.addCell("illumina Sys Error");
		table.addCell(getContStr(sm, smsolo, FilterResult.illuminaSpecific));
		table.addCell(getContStr(ti, tisolo, FilterResult.illuminaSpecific));
		
		table.addCell("only softClip support");
		table.addCell(getContStr(sm, smsolo, FilterResult.softClip));
		table.addCell(getContStr(ti, tisolo, FilterResult.softClip));

		table.addCell("Mutation at same Cycle");
		table.addCell(getContStr(sm, smsolo,
				FilterResult.MutationAtSameCycle));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.MutationAtSameCycle));

		table.addCell("Fisher test fail p > "
				+ (float) KarkinosProp.Fisher_Thres_For_SNV_Detection);
		table.addCell(getContStr(sm, smsolo, FilterResult.FisherTestFail));
		table.addCell(getContStr(ti, tisolo, FilterResult.FisherTestFail));


		table.addCell("SNV near Indel");
		table.addCell(getContStr(sm, smsolo,
				FilterResult.NEARINDEL));
		table.addCell(getContStr(ti, tisolo,
				FilterResult.NEARINDEL));
		return table;
	}

	private static void addchrCounts(Map<String, CounterA> map, String key) {
		CounterA ca = null;
		if(map.containsKey(key)){
			ca = map.get(key);
		}else{
			ca = new CounterA();
		}
		ca.inc();
		map.put(key, ca);
		
	}

	private static void setSeriespaint(JFreeChart jfc) {
		CategoryPlot cp = jfc.getCategoryPlot();
		try {
			CategoryItemRenderer cir = cp.getRenderer();
			cir.setSeriesPaint(0, Color.RED.darker());
			cir.setSeriesPaint(1, Color.BLUE.darker());
			cir.setSeriesPaint(2, Color.GREEN.darker());
			cir.setSeriesPaint(3, Color.GRAY.darker());
			cir.setSeriesPaint(4, Color.CYAN.darker());

		} catch (Exception ex) {
		}
		try {
			BarRenderer br = (BarRenderer) cp.getRenderer();
			br.setShadowVisible(false);
			br.setBarPainter(new StandardBarPainter());

		} catch (Exception ex) {
		}
	}

	private static String getContStr(Map<Integer, CounterA> map, int passFilter) {

		String s1 = "0";
		if (map.containsKey(passFilter)) {
			s1 = map.get(passFilter).getCnt() + "";
		}
		return s1;

	}

	private static String getContStr(Map<Integer, CounterA> map,
			Map<Integer, CounterA> map2, int passFilter) {

		String s1 = "0";
		String s2 = "0";
		if (map.containsKey(passFilter)) {
			s1 = map.get(passFilter).getCnt() + "";
		}
		if (map2.containsKey(passFilter)) {
			s2 = map2.get(passFilter).getCnt() + "";
		}

		return s1 + "(" + s2 + ")";

	}

	private static Object getTRTable(DataSet dataset, PeaksInfo pi) {

		AnalyseDist analyseDist = dataset.getAnalyseDist();
		try {

			Table table = new Table(6);
			// row1
			table.addCell("source");
			table.addCell("tumor purity");
			table.addCell("s.d.");
			table.addCell("#total SNP");
			table.addCell("correl");
			table.addCell("method");
			// row2
			table.addCell("n=1");
			table.addCell(format(analyseDist.getTumorratioFromLOH().getTumorratio()));
			table.addCell(format(analyseDist.getTumorratioFromLOH().getSd() * 0.01));
			table.addCell(format(analyseDist.getTumorratioFromLOH().getNumber()));
			table.addCell(format(analyseDist.getTumorratioFromLOH().getCorrel()));
			table.addCell(analyseDist.getTumorratioFromLOH().getModeStr());
			
			// row2
			table.addCell("n=3");
			table.addCell(format(analyseDist.getTumorratioFromGAIN().getTumorratio()));
			table.addCell(format(analyseDist.getTumorratioFromGAIN().getSd() * 0.01));
			table.addCell(format(analyseDist.getTumorratioFromGAIN().getNumber()));
			table.addCell(format(analyseDist.getTumorratioFromGAIN().getCorrel()));
			table.addCell(analyseDist.getTumorratioFromGAIN().getModeStr());
			
			// row2
			table.addCell("somatic");
			table.addCell(format(analyseDist.getTumorratioFromSomatic().getTumorratio()));
			table.addCell(format(analyseDist.getTumorratioFromSomatic().getSd() * 0.01));
			table.addCell(format(analyseDist.getTumorratioFromSomatic().getNumber()));
			table.addCell(format(analyseDist.getTumorratioFromSomatic().getCorrel()));
			table.addCell("from distribution of somatic mutations");
			
			table.addCell("tumor purity from ploidy matrix");
			table.addCell(format(dataset.getTumorratioFiitiingMatrix()));
			table.addCell("");
			table.addCell("");
			table.addCell("");
			table.addCell("");
			
			table.addCell("tumor purity from ploidy matrix");
			table.addCell(format(dataset.getTumorratioFiitiingMatrix()));
			table.addCell("");
			table.addCell("");
			table.addCell("");
			table.addCell("");
			
			table.addCell("ploidy");
			table.addCell(format(pi.getPloidy()));
			table.addCell("");
			table.addCell("");
			table.addCell("");
			table.addCell("");
			
			table.addCell("tumor purity used");
			table.addCell(format(dataset.getTumorRatio()));
			//
			
			
			return table;

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	private static JFreeChart getHetroSNPGraph(IntervalXYDataset normal,
			IntervalXYDataset tumor, float degree, float observedratio) {

		
		NumberAxis xAxis = new NumberAxis("ratio");
		xAxis.setRange(0, 100);
		NumberAxis yAxis = new NumberAxis("counts");
		yAxis.setAutoRangeIncludesZero(true);

		// make a horizontally combined plot
		// make a horizontally combined plot
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);

		final XYBarRenderer renderer0 = new XYBarRenderer();
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setShadowVisible(false);
		renderer0.setBarPainter(new StandardXYBarPainter());
		
		// add subplot 1...

		XYPlot subplot1 = new XYPlot(normal, xAxis, yAxis, renderer0);
		//subplot1.setDomainCrosshairVisible(true);
		//subplot1.setRangeCrosshairVisible(true);
		parent.add(subplot1, 1);
		setNoshawow(subplot1);

		// add subplot 2...
		final XYBarRenderer renderer1 = new XYBarRenderer();
		renderer1.setSeriesPaint(0, ChartColor.BLUE);
		renderer1.setShadowVisible(false);
		renderer1.setBarPainter(new StandardXYBarPainter());

		XYPlot subplot2 = new XYPlot(tumor, xAxis, yAxis, renderer1);
		//subplot2.setDomainCrosshairVisible(true);
		//subplot2.setRangeCrosshairVisible(true);
		parent.add(subplot2, 1);
		setNoshawow(subplot2);
		addMaker(subplot2, observedratio);
		addMaker(subplot2, (100 - observedratio));
		
		

		// createLineChart(data, string, string2, string3);
		int n = (int) (degree * 2);
		JFreeChart chart = new JFreeChart("Mutation allele frequency N=" + n,
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);
		
		chart.setBackgroundPaint(Color.WHITE);
		return chart;
	}

	private static void setNoshawow(XYPlot sp) {
		
		try {
			BarRenderer br = (BarRenderer) sp.getRenderer();
			br.setShadowVisible(false);
			br.setBarPainter(new StandardBarPainter());

		} catch (Exception ex) {
		}
		
	}

	private static void addMaker(XYPlot xyplot, double d) {
		Marker marker0 = new ValueMarker(d);
		marker0.setPaint(Color.RED);
		marker0.setOutlineStroke(new BasicStroke(5.0f));
		xyplot.addDomainMarker(marker0);

	}

	private static JFreeChart getHetroSNPGraph(IntervalXYDataset normal,
			IntervalXYDataset tumor, float degree) {

		
		NumberAxis xAxis = new NumberAxis("ratio");
		xAxis.setRange(0, 100);
		NumberAxis yAxis = new NumberAxis("counts");
		yAxis.setAutoRangeIncludesZero(true);

		// make a horizontally combined plot
		// make a horizontally combined plot
		CombinedRangeXYPlot parent = new CombinedRangeXYPlot(yAxis);

		final XYBarRenderer renderer0 = new XYBarRenderer();
		renderer0.setSeriesPaint(0, ChartColor.BLUE);
		renderer0.setShadowVisible(false);
		renderer0.setBarPainter(new StandardXYBarPainter());

		// add subplot 1...

		XYPlot subplot1 = new XYPlot(normal, xAxis, yAxis, renderer0);
		subplot1.setDomainCrosshairVisible(true);
		subplot1.setRangeCrosshairVisible(true);
		parent.add(subplot1, 1);

		// add subplot 2...
		final XYBarRenderer renderer1 = new XYBarRenderer();
		renderer1.setSeriesPaint(0, ChartColor.BLUE);
		renderer1.setShadowVisible(false);
		renderer1.setBarPainter(new StandardXYBarPainter());

		XYPlot subplot2 = new XYPlot(tumor, xAxis, yAxis, renderer1);
		subplot2.setDomainCrosshairVisible(true);
		subplot2.setRangeCrosshairVisible(true);
		parent.add(subplot2, 1);

		String n = format((degree * 2));
		// createLineChart(data, string, string2, string3);
		JFreeChart chart = new JFreeChart("Mutation allele frequency N=" + n,
				JFreeChart.DEFAULT_TITLE_FONT, parent, true);

		return chart;
	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private static String format(int num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}
}
