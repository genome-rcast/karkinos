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
package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class FunctionRegression implements java.io.Serializable {

	List<double[]> l;
	int minbaitlen = 0;
	int baitmargin = 30;
	int unit = 0;
	List<Point2D> gcMedianList = null;
	RegressionInfo regInfo = null;

	public static void main(String[] arg) {

		List<Point2D> list = new ArrayList<Point2D>();
		list.add(new Point2D.Double(25, 0.35));
		list.add(new Point2D.Double(30, 0.50));
		list.add(new Point2D.Double(35, 0.65));
		list.add(new Point2D.Double(40, 0.75));
		list.add(new Point2D.Double(45, 0.85));
		list.add(new Point2D.Double(50, 0.90));
		list.add(new Point2D.Double(55, 1));
		list.add(new Point2D.Double(60, 1.05));
		list.add(new Point2D.Double(65, 1.1));
		list.add(new Point2D.Double(70, 1.15));

		FunctionRegression fr = new FunctionRegression();
		fr.calcRegression(list);
		System.out.println(fr.a);
		System.out.println(fr.b);
		System.out.println(fr.c);
		for (int n = 25; n <= 75; n = n + 5) {
			System.out.println(fr.function(n));
		}

	}

	public FunctionRegression(List<double[]> l, int minbaitlen) {

		this.l = l;
		this.minbaitlen = minbaitlen;
		unit = ((minbaitlen / 10) * 10) + baitmargin;
		analyse();
		baseAdjust = getBaseAdjust(l);
		
		regInfo = new RegressionInfo();
		regInfo.setA(a);
		regInfo.setB(b);
		regInfo.setC(c);
		regInfo.setD(d);
		regInfo.setBaseAdjust(baseAdjust);
		regInfo.setFuncflg(funcflg);
		regInfo.setGcMedianList(gcMedianList);
		//for
		System.out.println("regression info");
		System.out.println("a="+a+"b="+b+"c="+c+"d="+d);
		for (int n = 25; n <= 75; n = n + 5) {
			System.out.println(n+"\t"+function(n));
		}
	}

	private float getBaseAdjust(List<double[]> l) {
		
		float startval = 0.4f;
		float endval = 1.6f;
		//decide bin
		int binidx = 0;
		int maxidx = 0;
		int maxdatacount = 0;
		for (float f = startval; f <= endval-0.4f; f = f + 0.1f) {
			
			int datacount = 0;
			float binstart = f;
			float binend = f+0.4f;
			for(double[] da : l){
				int x = (int) Math.round(da[0]);
				// gc%
				int y = (int) Math.round(da[1]);
				// value
				double z = da[2];
				if (z == 1) {
					continue;
				}
				double  val = getGCAdjustedVal(y,z);
				if(val==1)continue;
				if((val>=binstart)&&(val<=binend)){
					datacount++;
				}
				
			}	
			if(maxdatacount<datacount){
				maxdatacount=datacount;
				maxidx = binidx;
			}
			
			binidx++;
		}
		
		float binstartval = (float) (startval+(0.1*maxidx));
		float binendval = binstartval + 0.4f;
		
		float fmin = startval;
		float sdmin = 100;
		for (float f =  binstartval; f <= binendval; f = f + 0.01f) {
			
			SummaryStatistics ss = new SummaryStatistics();
			for(double[] da : l){
				int x = (int) Math.round(da[0]);
				// gc%
				int y = (int) Math.round(da[1]);
				// value
				double z = da[2];
				if (z == 1) {
					continue;
				}
				double  val = getGCAdjustedVal(y,z);
				if(val==1)continue;
				//
				if((val>binendval)||(val<binstartval)){
					continue;
				}
				//
				double diff = Math.abs(f - val);
				ss.addValue(diff);
				
			}			
			
			System.out.println(f+"\t"+ss.getStandardDeviation());
			if(sdmin>ss.getStandardDeviation()){
				
				sdmin=(float) ss.getStandardDeviation();
				fmin = f;
				
			}
			
		}
		//
		if((fmin==binstartval)||(fmin==binendval)){
			return 0f;
		}else{
			return fmin-1;
		}
	}

	public FunctionRegression() {
		// TODO Auto-generated constructor stub
	}

	public RegressionInfo getRegInfo() {
		return regInfo;
	}

	Map<String, SummaryStatistics> adjustmap = new HashMap<String, SummaryStatistics>();
	Map<Integer, DistMedian> map = new HashMap<Integer, DistMedian>();

	private void analyse() {
		//
		// SummaryStatistics total = new SummaryStatistics();
		for (double[] da : l) {
			// length
			int x = (int) Math.round(da[0]);
			// gc%
			int y = (int) Math.round(da[1]);
			// value
			double z = da[2];
			if (z == 1) {
				continue;
			}
			String key = getGrid(x, y);
			SummaryStatistics ss = adjustmap.get(key);
			if (ss == null) {
				ss = new SummaryStatistics();
			}
			ss.addValue(z);
			// total.addValue(z);
			adjustmap.put(key, ss);

			//
			DistMedian ssreg = map.get(y);
			if (ssreg == null) {
				ssreg = new DistMedian();
			}
			ssreg.addValue(z);
			map.put(y, ssreg);

		}
		//
		totalmean = 1.0f;

		//
		List<Point2D> list = new ArrayList<Point2D>();
		for (int n = 10; n < 90; n++) {
			//
			DistMedian ssreg = map.get(n);
			if (ssreg != null) {
				double x = ((double) n);
				double y = ssreg.getDistributionMedian();
				if(x>0&&y>0){
					System.out.println(x+"_"+y);
					list.add(new Point2D.Double(x,y));
				}
			}
		}
		//
		calcRegression(list);
		gcMedianList = list;

	}

	int gridnum = 20;

	private String getGrid(int x, int y) {

		int n = 0;
		int m = 0;

		// x length

		n = x / unit;
		if (n >= gridnum)
			n = gridnum;

		// y gc%
		int unitP = 100 / gridnum;// 5
		m = y / unitP;
		if (m >= gridnum)
			m = gridnum;

		return n + "-" + m;
	}

	double totalmean = 0;

	// y=ae^bx+c
	double zmin = 0;
	double a = 0, b = 0, c = 0, d= 0;
	double baseAdjust=0;

	int useXmin = 25;
	int useXMax = 75;
	int funcflg = 0;
	public static final int Poly2 = 1;
	public static final int Log = 2;
	public static final int Exp = 3;
	public static final int Poly3 = 4;

	private double[] calcRegression(List<Point2D> list) {

		//debug
		for (Point2D point : list) {
			System.out.println("X=\t"+point.getX()+"\t Y=\t"+point.getY());
		}
		
		list = trimList(list);
		
		//1.regression by 2nd degree polynomial
		int n = 0;
		double X = 0;
		double X2 = 0;
		double X3 = 0;
		double X4 = 0;
		
		//for poly3
		double X5 = 0;
		double X6 = 0;		

		double Y = 0;
		double XY = 0;
		double X2Y = 0;
		
		//for poly3
		double X3Y = 0;

		// LSM to find 2nd degree eq
		for (Point2D point : list) {

			//
			double x0 = point.getX();
			double y0 = point.getY();
			double x2 = Math.pow(x0, 2);
			double x3 = Math.pow(x0, 3);
			double x4 = x2 * x2;
			double x5 = x2 * x3;
			double x6 = x3 * x3;

			X += x0;
			X2 += x2;
			X3 += x3;
			X4 += x4;
			// poly3
			X5 += x5;
			X6 += x6;
			
			n++;
			
			Y += y0;
			XY += x0 * y0;
			X2Y += x2 * y0;
			
			// poly3
			X3Y += x3 * y0;
			

		}
		// resolve equation
		double[][] vals = { { n, X, X2 }, { X, X2, X3 }, { X2, X3, X4 } };
		double[] rhs = { Y, XY, X2Y };

		RealMatrix rm = new Array2DRowRealMatrix(vals);

		double a1 = 0;
		double b1 = 0;
		double c1 = 0;

		try {
			DecompositionSolver solver = new LUDecompositionImpl(rm).getSolver();
			RealVector b = new ArrayRealVector(rhs);
			RealVector x = solver.solve(b);
			a1 = x.getEntry(2);
			b1 = x.getEntry(1);
			c1 = x.getEntry(0);
		} catch (Exception ex) {
			//ex.printStackTrace();
		}
		double z = 0;
		
		double diffpow = 0;
		int n_2 = 0;
		for (Point2D point : list) {

			double x0 = point.getX();
			double y0 = point.getY();

			double y = a1*Math.pow(x0,2)+(b1*x0)+c1;
			double diff = y0 - y;
			diffpow = diffpow + Math.pow(diff,2);
			n_2++;
			z = z + pow(diff, 2);

		}
		double diffpoly =z;
		double BIC2 = n_2* Math.log((diffpow/n_2))+2*Math.log(n_2);
		
		//for poly3
		double[][] vals3 = { { n, X, X2, X3 }, { X, X2, X3,X4 }, 
				{ X2, X3, X4,X5 }, { X3, X4, X5,X6 } };
		double[] rhs3 = { Y, XY, X2Y, X3Y};

		RealMatrix rm3 = new Array2DRowRealMatrix(vals3);

		double a1_3 = 0;
		double b1_3 = 0;
		double c1_3 = 0;
		double d1_3 = 0;

		try {
			DecompositionSolver solver = new LUDecompositionImpl(rm3).getSolver();
			RealVector b = new ArrayRealVector(rhs3);
			RealVector x = solver.solve(b);
			a1_3 = x.getEntry(3);
			b1_3 = x.getEntry(2);
			c1_3 = x.getEntry(1);
			d1_3 = x.getEntry(0);
		} catch (Exception ex) {
			//ex.printStackTrace();
		}
		double z_3 = 0;
		double diffpow_3 = 0;
		int n_3 = 0;
		for (Point2D point : list) {

			double x0 = point.getX();
			double y0 = point.getY();

			double y = a1_3*Math.pow(x0,3)+(b1_3*(x0*x0))+(c1_3*x0)+d1_3;
			double diff = y0 - y;
			z_3 = z_3 + pow(diff, 2);

			diffpow_3 = diffpow_3 + Math.pow(diff,2);
			n_3++;
		}
		double diffpoly_3 =z_3;
		double BIC3 = n_3* Math.log((diffpow_3/n_3))+3*Math.log(n_3);
		
		//2 regression by logalistic
		double InX = 0;
		double InX2 = 0;
		Y = 0;
		double InXY = 0;

		for (Point2D point : list) {

			double x0 = point.getX();
			double y0 = point.getY();
			double lnX0 = Math.log(x0);

			//
			Y += y0;
			InXY += lnX0 * y0;
			InX += lnX0;
			InX2 += pow(lnX0, 2);

		}

		n = list.size();
		double ea = (n * InXY - InX * Y) / (n * InX2 - InX * InX);
		double eb = (InX2 * Y - InXY * InX) / (n * InX2 - InX * InX);
		//
		z = 0;
		for (Point2D point : list) {

			double x0 = point.getX();
			double y0 = point.getY();

			double y = ea * Math.log(x0) + eb;
			double diff = y0 - y;
			z = z + pow(diff, 2);

		}
		double alog = ea;
		double blog = eb;
		double zminlog = z;

		//3 regression by exponental
		for (float z0 = -0.15f; z0 < 0.15; z0 = z0 + 0.01f) {
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

				if ((y0 - z0) < 0) {
					break;
				}

				double val = (y0 - z0);
				if (Math.abs(val) < 0.0001) {
					boolean posi = true;
					if (val < 0) {
						posi = false;
					}
					val = 0.0001;
					if (!posi) {
						val = -1 * val;
					}

				}
				double lnY0 = Math.log(y0 - z0);
				InY += lnY0;
				XInY += x0 * lnY0;
				double x2 = pow(x0, 2);
				X += x0;
				X2 += x2;

			}

			if (n != size) {
				continue;
			}
			eb = (n * XInY - X * InY) / (n * X2 - X * X);
			ea = (X2 * InY - XInY * X) / (n * X2 - X * X);
			ea = Math.exp(ea);

			z = 0;
			for (Point2D point : list) {

				double x0 = point.getX();
				double y0 = point.getY();

				double y = ea * Math.exp(eb * x0) + z0;
				double diff = y0 - y;
				z = z + pow(diff, 2);

			}
			if (zmin == 0 || zmin > z) {
				zmin = z;
				a = ea;
				b = eb;
				c = z0;

			}

		}
		System.out.println("zpoly="+"\t"+ diffpoly+"zmineop="+zmin+"\t"+"zminlog="+zminlog);
		
		funcflg = getFuncFlg(diffpoly,diffpoly_3,BIC2,BIC3,zminlog,zmin);
		
		if(funcflg == Poly2){
			System.out.println("2nd degree polynomial");
			a = a1;
			b = b1;
			c= c1;
			return new double[] { a, b, c };
		}else if (funcflg == Poly3){
			
			a = a1_3;
			b = b1_3;
			c= c1_3;
			d = d1_3;
			return new double[] { a, b, c, d };
				
		}else if (funcflg == Exp) {
			System.out.println("exp");
			return new double[] { a, b, c };
		} else if(funcflg == Log){
			System.out.println("log");
			a = alog;
			b = blog;
			c=0;
			return new double[] { alog, blog ,0};
		}
		return null;

	}

//	public static final int Poly2 = 1;
//	public static final int Log = 2;
//	public static final int Exp = 3;
	private int getFuncFlg(double diffpoly,double diffpoly_3,
			double BIC2, double BIC3,
			double zminlog, double zminexp) {
		
		//
		if(diffpoly==0){
			diffpoly=10000;
		}
		if(diffpoly_3==0){
			diffpoly_3=10000;
		}
		if(zminlog==0){
			zminlog=10000;
		}
		if(zminexp==0){
			zminexp=10000;
		}
		if(BIC2<BIC3){
			if(diffpoly<zminlog && diffpoly<zminexp){
				return Poly2;
			}	
			if(zminlog<diffpoly && zminlog<zminexp){
				return Log;
			}	
			if(zminexp<zminlog && zminexp<diffpoly){
				return Exp;
			}	
			return Log;
		}else{
			if(diffpoly_3<zminlog && diffpoly_3<zminexp){
				return Poly3;
			}	
			if(zminlog<diffpoly_3 && zminlog<zminexp){
				return Log;
			}	
			if(zminexp<zminlog && zminexp<diffpoly_3){
				return Exp;
			}	
			return Log;
		}		
		
		
	}

	private List<Point2D> trimList(List<Point2D> list) {

		ArrayList al = new ArrayList<Point2D>();
		for (Point2D p : list) {
			if ((useXmin <= p.getX()) && (useXMax >= p.getX())) {
				al.add(p);
			}
		}
		return al;
	}

	private double pow(double x, int i) {
		return Math.pow(x, i);
	}

	private double getGCAdjustedVal(double y, double z) {

		double adjustmean = function(y);
		if(adjustmean<0.1){
			return z;
		}
		double diff = totalmean - adjustmean;
		z = z + diff;
		return z;

	}

	private double function(double x) {
		
		double val = 0;
		if (funcflg==Exp) {
			double y = a * Math.exp(b * x) + c;
			val =  y+baseAdjust;
		} else if(funcflg==Log) {
			double y = a * Math.log(x) + b;
			val = y+baseAdjust;
		}else if(funcflg==Poly2){
			double y = a * Math.pow(x, 2)+b*x+c;
			val = y+baseAdjust;
		}else if(funcflg==Poly3){
			double y = a * Math.pow(x, 3)+b*Math.pow(x, 2)+c*x+d;
			val = y+baseAdjust;
		}
		if(Double.isNaN(val)){
			return 0;
		}
		return val;
	}

	public double getAdjustedZ(double x, double y, double z) {

		if (z == 1)
			return 1;
		// GC% adjusted val
		z = getGCAdjustedVal(y, z);

		// remove high variance region

		// length
		int xi = (int) Math.round(x);
		// gc%
		int yi = (int) Math.round(y);
		String key = getGrid(xi, yi);
		SummaryStatistics ss = adjustmap.get(key);

		if (ss == null) {
			// no adjustment
			return trimBounds(z);
		}

//		if (ss.getN() >= 100 && ss.getStandardDeviation() > 2) {
//			// area from too high s.d
//			return step(z);			
//		}

		return trimBounds(z);
	}

//	private double step(double z) {
//
//		double a = Math.abs(z - 1.5);
//		double b = Math.abs(z - 1);
//		double c = Math.abs(z - 0.5);
//		if (a < b && a < c) {
//			return 1.5;
//		}
//
//		if (b < c && b < a) {
//			return 1;
//		}
//
//		if (c < b && c < a) {
//			return 1.5;
//		}
//
//		return z;
//	}

	private double trimBounds(double z) {
		if (z > 5) {
			z = 5;
		}
		if (z < 0.05) {
			z = 0.05;
		}
		return z;
	}

	public double[][] getMeanAry() {

		double[][] dary = new double[gridnum][gridnum];
		for (int n = 0; n < gridnum; n++) {
			for (int m = 0; m < gridnum; m++) {

				//
				String key = n + "-" + m;
				double val = 0;
				if (adjustmap.containsKey(key)) {
					val = adjustmap.get(key).getMean();
				}
				dary[n][m] = val;

			}
		}

		return dary;
	}

	public double[][] getSDAry() {

		double[][] dary = new double[gridnum][gridnum];
		for (int n = 0; n < gridnum; n++) {
			for (int m = 0; m < gridnum; m++) {

				//
				String key = n + "-" + m;
				double val = 0;
				if (adjustmap.containsKey(key)) {
					val = adjustmap.get(key).getStandardDeviation();
				}

				dary[n][m] = val;
			}
		}
		return dary;

	}

	public Object[] getCGLabel() {

		List<Integer> list = new ArrayList<Integer>();
		int dev = 100 / gridnum;
		for (int n = 0; n < gridnum; n++) {

			//
			list.add(dev * (n + 1));
		}
		return list.toArray();
	}

	public Object[] getBaitLabel() {
		List<Integer> list = new ArrayList<Integer>();
		int dev = 100 / gridnum;
		for (int n = 0; n < gridnum; n++) {

			//
			list.add(unit * (n + 1));
		}
		return list.toArray();
	}

}
