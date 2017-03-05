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
import java.util.Map;
import java.util.Properties;

import com.lowagie.text.Table;

import jp.ac.utokyo.rcast.karkinos.distribution.AnalyseDist;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.CounterA;
import jp.ac.utokyo.rcast.karkinos.readssummary.DepthCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsCounter;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;

public class CNVVcf {

	// ##fileformat=VCFv4.1
	// ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
	// ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
	// ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
	// ##ALT=<ID=CNV,Description="Copy number variable region">
	// ##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">

	public static void outData(String outpath, DataSet dataset) {

		//
		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			// header
			bw.write("#chr" + "\t");
			bw.write("start" + "\t");
			bw.write("end" + "\t");
			bw.write("copyNumber" + "\t");
			bw.write("gain/loss" + "\t");
			bw.write("num hetro SNP" + "\t");
			bw.write("AAF");
			bw.write("BAF");
			bw.write("type");
			
			bw.write("\n");


			writeCNV(bw, dataset, 1);
			//writeCNV(bw, dataset, 2);
			writeCNV(bw, dataset, 3);

			bw.close();
		} catch (IOException ex) {

		}

	}

	private static void writeCNV(BufferedWriter bw, DataSet dataset, int dispFlg) {

		try {

			
			for (CopyNumberInterval cni : dataset
					.getCopyNumberIntervalList(9)) {

//				if (cni.getCopynumber() == 2)
//					continue;
//				if (dispFlg == 1) {
//					if (cni.isAllelic())
//						continue;
////					if (cni.getCopynumber() == 0)
////						continue;
////					if (cni.getCopynumber() >= 4)
////						continue;
//				} else if (dispFlg == 2) {
//					if (!cni.isAllelic())
//						continue;
////					if (cni.getCopynumber() == 0)
////						continue;
////					if (cni.getCopynumber() >= 4)
////						continue;
//				} else {
//					if (1 <= cni.getCopynumber() && cni.getCopynumber() <= 4) {
//						continue;
//					}
//				}
				
				if (cni.getCopynumber() == dataset.getBaseploidy()){
					if(cni.getAaf()==dataset.getBaseploidy()/2){
						continue;
					}
				}
					
				if (dispFlg == 1) {
					if (cni.isAllelic())
						continue;
//					if (cni.getCopynumber() == 0)
//						continue;
//					if (cni.getCopynumber() >= 4)
//						continue;
				} else if (dispFlg == 2) {
					if (!cni.isAllelic())
						continue;
					if (cni.getCopynumber() == 0)
						continue;
					if (cni.isHdelation())
						continue;
				} else {
					if (1 <= cni.getCopynumber() && !cni.isHdelation()) {
						continue;
					}
					if(cni.getCopynumber()!=0 && cni.isHdelation()){
						continue;
					}
				}

				bw.write(cni.getChr() + "\t");
				bw.write(format(cni.getStart()) + "\t");
				bw.write(format(cni.getEnd()) + "\t");
				bw.write("n=" + cni.getCopynumber() + "\t");
				String s = cni.getCopynumber() < 2 ? "loss" : "gain";
				bw.write(s + "\t");
				
				bw.write(format(cni.getNoSNP()) + "\t");
				
				bw.write(format(cni.getAaf()) + "\t");
				bw.write(format(cni.getBaf()) + "\t");
				
				if (dispFlg == 1) {
					
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
