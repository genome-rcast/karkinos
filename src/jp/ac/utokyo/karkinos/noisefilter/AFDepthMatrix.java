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
package jp.ac.utokyo.karkinos.noisefilter;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;

public class AFDepthMatrix {

	int[][] counter = null;
	List<BinData> binlist;
	int binsize = 20;

	public AFDepthMatrix(List<BinData> binlist, int binsize) {

		this.binlist = binlist;
		this.binsize = binsize;

	}

	public static void main(String[] arg) {

		AFDepthMatrix inst = new AFDepthMatrix(null, 0);
		inst.a1 = 1.51;
		inst.b1 = -1.68;
		// double a1fix = 1.51;
		// double b1fix = -1.68;
		// for(double x=0;x<0.5;x=x+0.1){
		// System.out.println(x+"\t"+inst.func((float) x));
		// }
		float pval = inst.getPval(0.1f, 1);
		System.out.println(pval);

	}

	public float getPval(float adjusttedAF, float readdepth) {

		try {

			double x = Math.log(readdepth - c1) / (Math.log(a1) * b1);
			double sd = Math.abs(x * 0.5);
			if (x < 0) {
				return 0;
			}
			NormalDistribution nd = new NormalDistributionImpl(0, sd);
			return (float) (1 - nd.cumulativeProbability(adjusttedAF));

		} catch (Exception ex) {

		}
		return 1;
	}

	public void regSNV(SNVHolder snv, float tratio) {

		int depth = (int) (snv.getTumor().getTotalcnt() * tratio);
		double tr = CalcUtils.getTumorrateAdjustedRatio(snv, tratio);
		int depthindex = getDIndex(depth);
		//
		BinData bd = binlist.get(depthindex);
		//
		bd.setSNV(depth, tr);

	}

	public void sortList() {

		for (BinData bd : binlist) {
			bd.sortList();
		}

	}

	public void EMmethod(float tratio) {

		int idx = 0;
		for (BinData bd : binlist) {

			boolean firstBin = (idx == 0);
			boolean lastBin = (idx >= getMaxBinIdx(tratio));
			boolean lowdepthbin = (bd.getLdepth() < 20);

			if (firstBin || lastBin || lowdepthbin) {
				idx++;
				continue;
			}
			if (bd.getSamplesize() <= 10) {
				continue;
			}
			if (bd.getDepth() >= 100) {
				continue;
			}
			// excute EM method to separat noise to signal
			try {
				bd.em();
			} catch (Exception ex) {
				ex.printStackTrace();
			}
			idx++;
		}
		//
		SummaryStatistics ss = new SummaryStatistics();
		boolean atleast_one = false;
		for (BinData bd : binlist) {
			//
			if (bd.isEmExcuted()) {
				atleast_one = true;
				ss.addValue(bd.getNumCandidate());
			}
		}
		if (atleast_one) {
			for (BinData bd : binlist) {
				//
				bd.setPredictedCandnum((int) ss.getMean());
			}
		}
		//

	}

	private int getMaxBinIdx(float tratio) {
		if (binlist.size() == 0) {
			return 0;
		}
		BinData laste = binlist.get(binlist.size() - 1);
		int ld = laste.getLdepth();
		int dpethadj = Math.round((ld * tratio));
		int idx = 0;
		for (BinData bd : binlist) {

			//
			if (bd.includepeth(dpethadj)) {
				return idx;
			}
			idx++;
		}
		return idx;
	}

	public void regHetroSNP(SNVHolder snv, float tumorContentsRatio) {

		int depth = (int) snv.getTumor().getTotalcnt();
		double tr = getTR(snv, tumorContentsRatio);
		double secondAlleleF = getSecondAlleleF(snv, tumorContentsRatio);
		int depthindex = getDIndex(depth);

		BinData bd = binlist.get(depthindex);
		//
		bd.setHetroSNPAF(depth, tr);
		bd.setSecondAF(depth, secondAlleleF);

	}

	private double getTR(SNVHolder snv, float tumorContentsRatio) {

		double tAF = snv.getTumor().getRatio();
		double diffr = (0.5 - tAF) / tumorContentsRatio;
		return diffr + 0.5;

	}

	private double getSecondAlleleF(SNVHolder snv, float tumorContentsRatio) {

		double tAF = snv.getTumor().getSecondRatio();
		double diffr = (0.5 - tAF) / tumorContentsRatio;
		return diffr + 0.5;

	}

	private int getDIndex(int depth) {

		int ret = 0;
		for (int n = 0; n < binlist.size(); n++) {
			int deptht = binlist.get(n).getUdepth();
			if (depth < deptht) {
				ret = n;
				break;
			}
		}
		if (ret < 0)
			ret = 0;
		if (ret >= binsize)
			ret = binsize - 1;
		return ret;

	}

	public void calcregresion() {

		// border cond
		List<Point2D> raiodapthlist = new ArrayList<Point2D>();
		for (BinData bd : binlist) {
			//
			if (bd.getDepth() <= 100) {
				raiodapthlist.add(new Point2D.Double(bd.getAFBorder(), bd
						.getDepth()));
			}
		}
		calcRegression(raiodapthlist);
	}

	double a1 = 0;
	double b1 = 0;
	double c1 = 0;

	private void calcRegression(List<Point2D> list) {

		// 1.regression by 2nd degree polynomial
		// int n = 0;
		double X = 0;
		double X2 = 0;
		// double X3 = 0;
		// double X4 = 0;
		//
		//
		// double Y = 0;
		// double XY = 0;
		// double X2Y = 0;

		if (list.size() <= 3)
			return;

		// // LSM to find 2nd degree eq
		// for (Point2D point : list) {
		//
		// //
		// double x0 = point.getX();
		// double y0 = point.getY();
		// double x2 = Math.pow(x0, 2);
		// double x3 = Math.pow(x0, 3);
		// double x4 = x2 * x2;
		//
		// X += x0;
		// X2 += x2;
		// X3 += x3;
		// X4 += x4;
		//
		//
		// n++;
		//
		// Y += y0;
		// XY += x0 * y0;
		// X2Y += x2 * y0;
		//
		// }
		// // resolve equation
		//
		// double[][] vals = { { n, X, X2 }, { X, X2, X3 }, { X2, X3, X4 } };
		// double[] rhs = { Y, XY, X2Y };
		//
		// RealMatrix rm = new Array2DRowRealMatrix(vals);
		// try {
		// DecompositionSolver solver = new LUDecompositionImpl(rm).getSolver();
		// RealVector b = new ArrayRealVector(rhs);
		// RealVector x = solver.solve(b);
		// a1 = x.getEntry(2);
		// b1 = x.getEntry(1);
		// c1 = x.getEntry(0);
		// } catch (Exception ex) {
		// //ex.printStackTrace();
		// }
		int n = list.size();
		double zmin = Integer.MAX_VALUE;
		// 3 regression by exponental
		for (float z0 = 0f; z0 < 0.1; z0 = z0 + 0.01f) {
			//
			X = 0;
			X2 = 0;
			double InY = 0;
			double XInY = 0;

			int size = 0;
			for (Point2D point : list) {
				size++;
				double x0 = point.getX();
				double y0 = point.getY();
				x0 = x0 - z0;

				double lnY0 = Math.log(y0 - z0);
				InY += lnY0;
				XInY += x0 * lnY0;
				double x2 = Math.pow(x0, 2);
				X += x0;
				X2 += x2;

			}

			if (n != size) {
				continue;
			}
			double eb = (n * XInY - X * InY) / (n * X2 - X * X);
			double ea = (X2 * InY - XInY * X) / (n * X2 - X * X);
			ea = Math.exp(ea);

			double z = 0;
			for (Point2D point : list) {

				double x0 = point.getX();
				double y0 = point.getY();

				double y = ea * Math.exp(eb * (x0 - z0));

				double diff = y0 - y;
				z = z + Math.pow(diff, 2);

			}
			if (zmin == 0 || zmin > z) {
				zmin = z;
				a1 = ea;
				b1 = eb;
				c1 = z0;

			}

		}

	}

	public double func(float x) {

		// double y = _func(x);
		// double localmin = (-1*b1)/(2*a1);
		// if(x>localmin){
		// y = _func((float)localmin);
		// }
		if ((a1 == 0) && (b1 == 0) && (c1 == 0) || (b1 > 0)) {
			double y = a1fix * Math.exp(b1fix * x);
			if (y < 0)
				y = 0;
			return y;
		}

		double y = a1 * Math.exp(b1 * (x - c1));
		if (y < 0)
			y = 0;
		return y;
	}

	private float ratio(int total2, int total) {

		return (float) ((double) (total2) / (double) (total));
	}

	private void printDist() {

		for (int n = 0; n < 20; n++) {

			for (int m = 0; m < 20; m++) {
				System.out.print(counter[m][n] + "\t");
			}
			System.out.println();
		}

	}

	public boolean reject(int depth_c, float adjusttedAF,
			boolean highErrorSample) {

		boolean reject = _reject(depth_c, adjusttedAF);
		
		
		if (!reject) {

			if ((depth_c <= 100)
					&& (adjusttedAF < KarkinosProp.mintumorratioForFilter1)) {
				reject = true;
			}

		}
		if (highErrorSample) {
			if (reject) {
				if (adjusttedAF > 0.6 && depth_c >= 30) {
					reject = false;
				}
				if (adjusttedAF > 0.8) {
					reject = false;
				}
			}
		} else {
			if (reject) {
				if (adjusttedAF > 0.2 && depth_c >= 10) {
					reject = false;
				}
				if (adjusttedAF > 0.3) {
					reject = false;
				}
			}
		}
		return reject;

	}

	public boolean _reject(int depth_c, float adjusttedAF) {

		if ((a1 == 0) && (b1 == 0) && (c1 == 0) || (b1 > 0)) {
			return defultreject(depth_c, adjusttedAF);
		}

		double y = func(adjusttedAF);
		return depth_c < y;

		// //case 1 no list
		// if((raiodapthlist==null)||(raiodapthlist.size()==0)){
		// return defultreject(depth_c,adjusttedAF);
		// }
		// double x = getTargetX(depth_c,raiodapthlist);
		// if(adjusttedAF<=x){
		// return true;
		// }else{
		// return false;
		// }

	}

	static final double a1fix = 1.51;
	static final double b1fix = -1.68;

	public static boolean defultreject(int depth_c, float x) {

		// y=-500x+200
		// double y = (-500*adjusttedAF)+200;
		// if(y<0)y=0;
		// if(depth_c<=20){
		// return true;
		// }

		double y = a1fix * Math.exp(b1fix * x);
		if (y < 0)
			y = 0;
		return depth_c < y;

	}

	private double getTargetX(int depth_c, List<Point2D> list) {

		Point2D laste = list.get(list.size() - 1);
		double lastd = laste.getY();
		if (depth_c > lastd) {
			return laste.getX();
		}
		if (list.size() == 1) {
			return list.get(0).getX();
		}

		for (int n = 0; n < list.size() - 1; n++) {

			Point2D e1 = list.get(n);
			Point2D e2 = list.get(n + 1);
			//
			if ((depth_c >= e1.getY()) && (depth_c <= e2.getY())) {
				//
				double mean = (e1.getX() + e2.getX()) / 2;
				return mean;
			}
		}
		return list.get(0).getX();

	}

}
