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
package jp.ac.utokyo.rcast.karkinos.graph.output;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Map.Entry;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.DepthCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;

public class TextSummary {
	public static void outTextData(String outpath, ReadsSummary readsSummary,
			String readsStat, DataSet dataset, float ploidy) {
		try {
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			writeReadSummary(bw, readsSummary, readsStat);
			writeTR(bw, dataset, ploidy);
			writeSNPCorrel(bw, dataset);
			writeCNV(bw, dataset, 1);
			writeCNV(bw, dataset, 2);
			writeCNV(bw, dataset, 3);

			bw.close();
		} catch (IOException ex) {
		}
	}

	private static void writeSNPCorrel(BufferedWriter bw, DataSet dataset)
			throws IOException {
		int snpRecurrent = 0;
		int snpDiff = 0;

		List<Double> xl = new ArrayList<Double>();
		List<Double> yl = new ArrayList<Double>();

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

							xl.add(x);
							yl.add(y);

							if ((x >= 0.4)) {
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

		bw.write("recurrent SNP=\t" + snpRecurrent + "\n");
		bw.write("different SNP=\t" + snpDiff + "\n");
		double parcentDifferent = (double) snpDiff
				/ (double) (snpRecurrent + snpDiff);
		bw.write("diff SNP %=\t" + parcentDifferent + "\n");

		double correl = correlation(xl, yl);
		bw.write("correl=\t" + correl + "\n");
	}

	public static double correlation(List<Double> xs, List<Double> ys) {
		// TODO: check here that arrays are not null, of the same length etc

		double sx = 0.0;
		double sy = 0.0;
		double sxx = 0.0;
		double syy = 0.0;
		double sxy = 0.0;

		int n = xs.size();

		for (int i = 0; i < n; ++i) {
			double x = xs.get(i);
			double y = ys.get(i);

			sx += x;
			sy += y;
			sxx += x * x;
			syy += y * y;
			sxy += x * y;
		}

		// covariation
		double cov = sxy / n - sx * sy / n / n;
		// standard error of x
		double sigmax = Math.sqrt(sxx / n - sx * sx / n / n);
		// standard error of y
		double sigmay = Math.sqrt(syy / n - sy * sy / n / n);

		// correlation is just a normalized covariation
		return cov / sigmax / sigmay;
	}

	private static void writeReadSummary(BufferedWriter bw,
			ReadsSummary readsSummary, String readsStat) {
		try {
			ReadsCounter normal = readsSummary.getNormalCounter();
			ReadsCounter tumor = readsSummary.getTumorCounter();

			// row1
			bw.write("#reads stats" + "\t");
			bw.write("normal" + "\t");
			bw.write("tumor" + "\t");
			bw.write("\n");

			// // row
			if (readsStat != null) {
				String[] sa = readsStat.split(",");
				String normalfile = sa[0];
				Properties np = getProp(normalfile);
				String tumorfile = sa[1];
				Properties tp = getProp(tumorfile);

				String[] keys = new String[] { "total", "passfilter", "mapped",
						"nonduplicate", "duplicate", "unique", "proper",
						"properOrunique", "afterrealgin", "identityLow" };
				for (String key : keys) {
					if (np.containsKey(key)) {
						// row2
						String kyeDisp = key.replaceAll("Or", " or ");
						bw.write(kyeDisp + " tags" + "\t");
						bw.write(format(np.getProperty(key)) + "\t");
						bw.write(format(tp.getProperty(key)) + "\t");
						bw.write("\n");

						if (key.equals("duplicate")) {
							try {
								// row 13
								bw.write("duplicate (%) \t");
								bw.write(format(getParcent(np.getProperty(key),
										np.getProperty("nonduplicate")))
										+ "%"
										+ "\t");
								bw.write(format(getParcent(tp.getProperty(key),
										tp.getProperty("nonduplicate")))
										+ "%"
										+ "\t");
								bw.write("\n");
							} catch (Exception ex) {
							}
						}
					}
				}
			}

			// row
			bw.write("Total used tags" + "\t");
			bw.write(format(normal.getTotalmap()) + "\t");
			bw.write(format(tumor.getTotalmap()) + "\t");
			bw.write("\n");

			// row4
			bw.write("Target tags" + "\t");
			bw.write(format(normal.getTotalOnTarget()) + "\t");
			bw.write(format(tumor.getTotalOnTarget()) + "\t");
			bw.write("\n");
			// row5
			bw.write("Target unique tags" + "\t");
			bw.write(format(normal.getTotalUniqueOntarget()) + "\t");
			bw.write(format(tumor.getTotalUniqueOntarget()) + "\t");
			bw.write("\n");

			if (readsSummary.isPairStats()) {
				// row6
				bw.write("First reads" + "\t");
				bw.write(format(normal.getFirstReads()) + "\t");
				bw.write(format(tumor.getFirstReads()) + "\t");
				bw.write("\n");
				// row7
				bw.write("Second reads" + "\t");
				bw.write(format(normal.getSecondReads()) + "\t");
				bw.write(format(tumor.getSecondReads()) + "\t");
				bw.write("\n");
				// row8
				bw.write("Both mapped" + "\t");
				bw.write(format(normal.getBothmap()) + "\t");
				bw.write(format(tumor.getBothmap()) + "\t");
				bw.write("\n");
				// row9
				bw.write("Proper reads" + "\t");
				bw.write(format(normal.getPropermap()) + "\t");
				bw.write(format(tumor.getPropermap()) + "\t");
				bw.write("\n");
				// row10
				bw.write("mean library size" + "\t");
				bw.write(format(normal.getMeanInsertSize()) + "\t");
				bw.write(format(tumor.getMeanInsertSize()) + "\t");
				bw.write("\n");
			}

			// // row11
			// row 12
			bw.write("mean depth extended" + "\t");
			bw.write(format(readsSummary.getNormalDepth().getMeanDepth())
					+ "\t");
			bw.write(format(readsSummary.getTumorDepth().getMeanDepth()) + "\t");
			bw.write("\n");

			bw.write("mean depth" + "\t");
			bw.write(format(readsSummary.getNormalDepth().getMeanOntagDepth())
					+ "\t");
			bw.write(format(readsSummary.getTumorDepth().getMeanOntagDepth())
					+ "\t");
			bw.write("\n");
			// row 12
			bw.write("On target(%)" + "\t");
			bw.write(format(normal.onTargetParcent()) + "%" + "\t");
			bw.write(format(tumor.onTargetParcent()) + "%" + "\t");
			bw.write("\n");

			bw.write("more than X20(%)" + "\t");
			bw.write(format(getX20percent(readsSummary.getNormalDepth())) + "%"
					+ "\t");
			bw.write(format(getX20percent(readsSummary.getTumorDepth())) + "%"
					+ "\t");
			bw.write("\n");

			bw.write("mean depth CDS" + "\t");
			bw.write(format(readsSummary.getNormalDepth().getMeanCDSDepth())
					+ "\t");
			bw.write(format(readsSummary.getTumorDepth().getMeanCDSDepth())
					+ "\t");
			bw.write("\n");

			bw.write("more than X10(%) CDS" + "\t");
			bw.write(format(readsSummary.getNormalDepth().getOver10X()) + "%"
					+ "\t");
			bw.write(format(readsSummary.getTumorDepth().getOver10X()) + "%"
					+ "\t");
			bw.write("\n");

			bw.write("more than X20(%) CDS" + "\t");
			bw.write(format(readsSummary.getNormalDepth().getOver20X()) + "%"
					+ "\t");
			bw.write(format(readsSummary.getTumorDepth().getOver20X()) + "%"
					+ "\t");
			bw.write("\n");

			Map<String, ReadsCounter> datamapn = readsSummary
					.getNormalPerChrom();
			Map<String, ReadsCounter> datamapt = readsSummary
					.getTumorPerChrom();

			bw.write("#reads counts by chr" + "\n");
			for (Entry<String, ReadsCounter> entry : datamapt.entrySet()) {
				String chrom = entry.getKey();
				bw.write("ReadsCounts by " + chrom + "\t");
				ReadsCounter rct = entry.getValue();
				ReadsCounter rcn = datamapn.get(chrom);
				int nmap = 0;
				if (rcn != null) {
					nmap = rcn.getTotalmap();
				}
				bw.write(format(rct.getTotalmap()) + "\t");
				bw.write(format(nmap) + "\t");
				bw.write("\n");
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static double getParcent(String dup, String nondup) {
		double d = Double.parseDouble(dup);
		double nond = Double.parseDouble(nondup);
		d = (d / (d + nond)) * 100;
		return d;
	}

	private static void writeTR(BufferedWriter bw, DataSet dataset, float ploidy) {
		try {
			AnalyseDist dist = dataset.getAnalyseDist();
			// row1
			bw.write("#source" + "\t");
			bw.write("tumor_rate" + "\t");
			bw.write("s.d." + "\t");
			bw.write("#total SNP" + "\t");
			bw.write("correl" + "\t");
			bw.write("method" + "\t");
			bw.write("\n");
			// row2
			bw.write("n=1" + "\t");
			bw.write(format(dist.getTumorratioFromLOH().getTumorratio()) + "\t");
			bw.write(format(dist.getTumorratioFromLOH().getSd() * 0.01) + "\t");
			bw.write(format(dist.getTumorratioFromLOH().getNumber()) + "\t");
			bw.write(format(dist.getTumorratioFromLOH().getCorrel()) + "\t");
			bw.write(dist.getTumorratioFromLOH().getModeStr() + "\t");
			bw.write("\n");
			// row2
			bw.write("n=3" + "\t");
			bw.write(format(dist.getTumorratioFromGAIN().getTumorratio())
					+ "\t");
			bw.write(format(dist.getTumorratioFromGAIN().getSd() * 0.01) + "\t");
			bw.write(format(dist.getTumorratioFromGAIN().getNumber()) + "\t");
			bw.write(format(dist.getTumorratioFromGAIN().getCorrel()) + "\t");
			bw.write(dist.getTumorratioFromGAIN().getModeStr() + "\t");
			bw.write("\n");

			// somatic
			bw.write("somatic mutation" + "\t");
			bw.write(format(dist.getTumorratioFromSomatic().getTumorratio())
					+ "\t");
			bw.write(format(dist.getTumorratioFromSomatic().getSd() * 0.01)
					+ "\t");
			bw.write(format(dist.getTumorratioFromSomatic().getNumber()) + "\t");
			// bw.write(format(dist.getTumorratioFromSomatic().getCorrel())+"\t");
			// bw.write(dist.getTumorratioFromSomatic().getModeStr()+"\t");
			bw.write("\n");

			//
			bw.write("tc used" + "\t");
			bw.write(format(dataset.getTumorRatio()) + "\t");
			bw.write(format(dist.getTcflg()) + "\t");
			bw.write("\n");

			bw.write("ploidy" + "\t");
			bw.write(format(ploidy) + "\t");
			bw.write(format(dist.getTcflg()) + "\t");
			bw.write("\n");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void writeCNV(BufferedWriter bw, DataSet dataset, int dispFlg) {
		try {
			// header
			bw.write("#chr" + "\t");
			bw.write("start" + "\t");
			bw.write("end" + "\t");
			bw.write("copy number" + "\t");
			bw.write("gain/loss" + "\t");
			if (dispFlg == 1) {
				bw.write("num hetro SNP" + "\t");
				bw.write("AAF" + "\t");
				bw.write("BAF" + "\t");
			}
			bw.write("type");

			bw.write("\n");

			for (CopyNumberInterval cni : dataset.getCopyNumberIntervalList(9)) {
				if (cni.getCopynumber() == dataset.getBaseploidy()) {
					if (cni.getAaf() == dataset.getBaseploidy() / 2) {
						continue;
					}
				}
				if (dispFlg == 1) {
					if (cni.isAllelic())
						continue;
					// if(cni.getCopynumber()==0)continue;
					// if(cni.getCopynumber()>=4)continue;
				} else if (dispFlg == 2) {
					if (!cni.isAllelic())
						continue;
					if (cni.getCopynumber() == 0)
						continue;
					if (cni.isHdeletion())
						continue;
				} else {
					if (1 <= cni.getCopynumber() && !cni.isHdeletion()) {
						continue;
					}
					if (0.5 <= cni.getCopynumber() && cni.isHdeletion()) {
						continue;
					}
					if (0 > cni.getCopynumber() && cni.isHdeletion()) {
						continue;
					}
				}

				bw.write(cni.getChr() + "\t");
				bw.write(format(cni.getStart()) + "\t");
				bw.write(format(cni.getEnd()) + "\t");
				bw.write("n=" + cni.getCopynumber() + "\t");
				String s = cni.getCopynumber() < 2 ? "loss" : "gain";
				bw.write(s + "\t");
				if (dispFlg == 1) {
					bw.write(format(cni.getNoSNP()) + "\t");
					// bw.write(cni.isSupportbyAllelic()+"\t");
					bw.write(cni.getAaf() + "\t");
					bw.write(cni.getBaf() + "\t");
					bw.write("fromDepth");
				} else if (dispFlg == 2) {
					bw.write("allelic");
				} else {
					if (cni.getCopynumber() < 1) {
						bw.write("HD");
					} else {
						bw.write("Amp");
					}
				}

				bw.write("\n");
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
			return (float) (((double) (total - less20x) / (double) total) * 100);
		} catch (Exception ex) {
		}
		return 0f;
	}
}
