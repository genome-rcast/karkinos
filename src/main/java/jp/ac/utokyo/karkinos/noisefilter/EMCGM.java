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

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class EMCGM {

	public static void main(String[] arg) {

		EMCGM emcgm = new EMCGM();
		emcgm.setInitVariance(0.05);
		List<Point2D> lista = new ArrayList<Point2D>();
		//
		lista.add(new Point2D.Double(0.01, 10));
		lista.add(new Point2D.Double(0.02, 10));
		lista.add(new Point2D.Double(0.1, 10));
		lista.add(new Point2D.Double(0.12, 10));
		lista.add(new Point2D.Double(0.101, 10));
		lista.add(new Point2D.Double(0.2, 10));
		lista.add(new Point2D.Double(0.25, 10));
		lista.add(new Point2D.Double(0.33, 10));
		lista.add(new Point2D.Double(0.39, 10));
		lista.add(new Point2D.Double(0.35, 10));
		lista.add(new Point2D.Double(0.4, 10));
		lista.add(new Point2D.Double(0.45, 10));
		lista.add(new Point2D.Double(0.41, 10));
		lista.add(new Point2D.Double(0.47, 10));
		lista.add(new Point2D.Double(0.5, 10));

		emcgm.analyse(lista);
	}

	public static final double u1Min = 0.35;
	public static final double u1Max = 0.75;

	public List<EMCGMBean> analyse(List<Point2D> list) {

		if(list==null||list.size()==0){
			return null;
		}
		int initnoisecnt = 0;
		for (Point2D p2d : list) {
			//System.out.println(p2d.getX());
			if(p2d.getX() < 0.3){
				initnoisecnt++;
			}
		}
		double initr = 0.5;
		if(initnoisecnt>0){
			initr = (double)initnoisecnt/(double)list.size();
		}
		double r = initr, v0 = initVariance, u0 = 0, v1 = initVariance, u1 = 0.5;

		int cntflg = 0;
		int maxiterate = 3;
		int cnt=0;
		List<EMCGMBean> dlist=null;
		while (cntflg == 0) {
			
			// estep
			dlist = estep(list, r, u0, v0, u1, v1);
			// mstep
			double[] d = mstep(dlist);
			cntflg = checkreality(d);
			if(cntflg!=0){
				break;
			}			
			r = d[0];
			u0 =d[1];
			v0 =d[2];
			u1 =d[3];
			v1 =d[4];		
			
			
			cnt++;
			if(cnt>=maxiterate){
				break;
			}
			if(NaN(d)){
				break;
			}
		}
		int candcnt = 0;
//		NormalDistribution noiseP = new NormalDistributionImpl(u0,Math.sqrt(v0));
		double borderAFmin = 0;
		double borderAFmax = 0;
//		try {
//			borderAF = noiseP.inverseCumulativeProbability(0.95);
//		} catch (MathException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}	
		for(EMCGMBean bean:dlist){
			
			//System.out.println(bean.getValue() +"\t" + bean.getP0());
			double psignalval = bean.getP1()/(bean.getP0()+bean.getP1());
			boolean signal = (psignalval>0.5); 
			if(signal){
				
				candcnt++;
				if((borderAFmax ==0) || bean.getValue()<borderAFmax){
					borderAFmax = bean.getValue();
				}				
				
			}else{
				
				if((borderAFmin ==0) || bean.getValue()>borderAFmin){
					borderAFmin = bean.getValue();
				}
				
			}
			System.out.println(bean.getValue() +"\t" + bean.getP0()+"\t" + borderAF +"\t"+signal);
			
		}
				
		double border = solveBorder(r,u0,v0,v1,u1);
		if(border<0){
			border = borderAFmin;
		}
		this.borderAF = border;
		numcandidate = candcnt;
		
		success = true;
		return dlist;
		
	}

	private double solveBorder(double r, double u0, double v0, double v1,
			double u1) {
		
		NDist noisepeaks = new NDist(u0, v0);
		NDist signalpeaks = new NDist(u1, v1);
		
		for (double x = u0;x<u1;x=x+0.01) {

			double value = x;
			double like0 = r * noisepeaks.getP(value);
			if(Double.isNaN(like0))like0=0;
			double like1 = (1 - r) * signalpeaks.getP(value);
			if(Double.isNaN(like1))like1=0;
			
			double p0 = like0 / (like0 + like1);
			double p1 = like1 / (like0 + like1);
			double pnoiseval = p0/(p0+p1);
			if(pnoiseval<0.95){
				return x;
			}
			
		}
		return -1;

	}

	private int checkreality(double[] d) {
		
		double u0 =d[1];
		double v0 =d[2];
		double u1 =d[3];
		double v1 =d[4];
		
		boolean v0toolow = (v0 <(initVariance/3));
		boolean v1toolow = (v1 <(initVariance/3));
		boolean v0toohigh = (v0 >(initVariance*3));
		boolean v1toohigh = (v1 >(initVariance*3));
		
		boolean u1toolow = u1 < u1Min;
		boolean u1toohigh = u1 < u1Max;
		
		if(v0toolow || v1toolow || v0toohigh || v1toohigh || u1toolow || u1toohigh){
			return 1;
		}		
		return 0;
		
	}

	private boolean NaN(double[] d) {
		
		for(double dv:d){
			if(Double.isNaN(dv)||Double.isInfinite(dv)){
				return true;
			}
		}
		return false;
	}

	private double[] mstep(List<EMCGMBean> dlist) {

		double[] d = new double[6];
		double r = 0, u0 = 0, v0 = 0, u1 = 0, v1 = 0;
		int cntflg = 0;
		//
		int n = 0;
		double sumZ0 = 0;
		double sumZ1 = 0;
		double sumZ0X = 0;
		double sumZ1X = 0;
		double sumZ0X2 = 0;

		for (EMCGMBean bean : dlist) {

			//
			n++;
			sumZ0 = sumZ0 + bean.getP0();
			sumZ1 = sumZ1 + bean.getP1();

			sumZ0X = sumZ0X + (bean.getP0() * bean.getValue());
			sumZ1X = sumZ1X + (bean.getP1() * bean.getValue());

			sumZ0X2 = sumZ0X2 + (bean.getP0() * pow2(bean.getValue()));

		}
		r = sumZ0 / (double) n;
		u0 = 0;// do not change
		u1 = sumZ1X / sumZ1;

		

		v0 = sumZ0X2 / sumZ0;

		double sumZ1X2 = 0;
		for (EMCGMBean bean : dlist) {

			//
			sumZ1X2 = sumZ1X2 + (bean.getP0() * pow2(bean.getValue() - u1));

		}
		v1 = sumZ1X2 / sumZ1;

		//
		d[0] = r;
		d[1] = u0;
		d[2] = v0;
		d[3] = u1;
		d[4] = v1;
		d[5] = cntflg;
		print(d);
		return d;

	}

	private void print(double[] d) {

		System.out.println(d[0] + "\t" + d[1] + "\t" + d[2] + "\t" + d[3]
				+ "\t" + d[4] + "\t" + d[5] + "\t");

	}

	private double pow2(double x) {
		return Math.pow(x, 2);
	}

	private List<EMCGMBean> estep(List<Point2D> list, double r, double u0,
			double v0, double u1, double v1) {

		// restrict u range
		if (u1 <= u1Min)
			u1 = u1Min;

		NDist noisepeaks = new NDist(u0, v0);
		NDist signalpeaks = new NDist(u1, v1);

		List<EMCGMBean> data = new ArrayList<EMCGMBean>();

		for (Point2D p2d : list) {

			double value = p2d.getX();
			double like0 = r * noisepeaks.getP(value);
			double like1 = (1 - r) * signalpeaks.getP(value);
			double p0 = like0 / (like0 + like1);
			double p1 = like1 / (like0 + like1);
			//
			EMCGMBean emcb = new EMCGMBean();
			emcb.setLike0(like0);
			emcb.setLike1(like1);
			emcb.setP0(p0);
			emcb.setP1(p1);
			emcb.setValue(value);
			data.add(emcb);
		}
		return data;

	}

	int numcandidate;
	int numnoise;
	double initVariance = 0.1;
	boolean success = false;

	public boolean isSuccess() {
		return success;
	}

	public int getNumcandidate() {
		return numcandidate;
	}

	public int getNumnoise() {
		return numnoise;
	}

	public double getBorderAF() {
		return borderAF;
	}

	public void setInitVariance(double initVariance) {
		this.initVariance = initVariance;
	}

	double borderAF;

	public void _analyse(List<Point2D> list) {

		//

	}

}
