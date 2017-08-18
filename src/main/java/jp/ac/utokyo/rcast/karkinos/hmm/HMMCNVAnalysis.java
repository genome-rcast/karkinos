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
package jp.ac.utokyo.rcast.karkinos.hmm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfIntegerFactory;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class HMMCNVAnalysis {

	public static void main(String[] aa) {

//		OpdfIntegerFactory factory = new OpdfIntegerFactory(4);
//		Hmm<ObservationInteger> hmm = new Hmm<ObservationInteger>(4, factory);
//
//		hmm.setPi(0, 0.05);
//		hmm.setPi(1, 0.90);
//		hmm.setPi(2, 0.04);
//		hmm.setPi(3, 0.01);
//
//		hmm.setOpdf(0, new OpdfInteger(new double[] { 0.9, 0.1, 0, 0 }));
//		hmm.setOpdf(1, new OpdfInteger(new double[] { 0.05, 0.9, 0.05, 0 }));
//		hmm.setOpdf(2, new OpdfInteger(new double[] { 0, 0.08, 0.9, 0.02 }));
//		hmm.setOpdf(3, new OpdfInteger(new double[] { 0, 0, 0.1, 0.9 }));
//
//		hmm.setAij(0, 1, 0.1);
//		hmm.setAij(0, 0, 0.9);
//		hmm.setAij(0, 2, 0);
//		hmm.setAij(0, 3, 0);
//
//		hmm.setAij(1, 0, 0.05);
//		hmm.setAij(1, 1, 0.9);
//		hmm.setAij(1, 2, 0.05);
//		hmm.setAij(1, 3, 0);
//
//		hmm.setAij(2, 0, 0);
//		hmm.setAij(2, 1, 0.08);
//		hmm.setAij(2, 2, 0.9);
//		hmm.setAij(2, 3, 0.2);
//
//		hmm.setAij(3, 0, 0);
//		hmm.setAij(3, 1, 0);
//		hmm.setAij(3, 2, 0.1);
//		hmm.setAij(3, 3, 0.9);
//
//		List<ObservationInteger> l = new ArrayList<ObservationInteger>();
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(0));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(0));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(2));
//		l.add(new ObservationInteger(2));
//		l.add(new ObservationInteger(0));
//		l.add(new ObservationInteger(2));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//		l.add(new ObservationInteger(1));
//
//		int[] ary = hmm.mostLikelyStateSequence(l);
//		for (int n : ary) {
//			System.out.println(n);
//		}

	}

	public static void calc(DataSet dataset) throws IOException {

		Hmm<ObservationReal> hmm = getHMM(dataset);
		List<List<WaveletIF>> plist = dataset.getCapInterval();
		for (List<WaveletIF> list : plist) {

			//
			List<ObservationReal> olist = getList(list);
			int[] hmmary = hmm.mostLikelyStateSequence(olist);
			int cnt = 0;
			for (int m : hmmary) {

				//
				WaveletIF wif = list.get(cnt);
				wif.setHMMValue((m + 1) * 0.5);
				cnt++;

			}

		}

	}

	private static List<ObservationReal> getList(List<WaveletIF> list) {

		List<ObservationReal> olist = new ArrayList<ObservationReal>();
		for (WaveletIF wi : list) {

			ObservationReal oi = new ObservationReal(wi.getDenioseValue());
			olist.add(oi);
		}
		return olist;
	}

	private static Hmm<ObservationReal> getHMM(DataSet dataset) {

		int[] countn = new int[10];
		SummaryStatistics ss = new SummaryStatistics();
		Set<Double> checkReg = new HashSet<Double>();
		List<List<WaveletIF>> plist = dataset.getCapInterval();
		for (List<WaveletIF> list : plist) {

			for (WaveletIF wi : list) {

				double cn = wi.getCN();
				int n = (int) (cn / 0.5);
				if (n >= 4)
					n = 4;
				countn[n - 1] = countn[n - 1] + 1;
				
				double dVal = wi.getDenioseValue();				
				if(!checkReg.contains(dVal)){
						//filter alreadey add val to calculate realistic variance
						//in wavelet denoized value list
						ss.addValue(dVal);
						checkReg.add(dVal);					
				}

			}

		}
		int nodesize = 0;
		int total = 0;
		for (int m : countn) {
			if (m > 0) {
				nodesize++;
			}
			total = total + m;
		}
		//asuume normal distribution for each n=1,2,3,4
		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(4, factory);
		int idx = 0;
		for (int m : countn) {

			double p = (double) m / (double) total;
			System.out.println(idx + "\t" + p);
			if (idx >= nodesize)
				break;
			hmm.setPi(idx, p);
			idx++;
		}
		int multifuctor = 2;
		double variance = ss.getVariance();
		variance = variance*multifuctor;
		System.out.println("variance="+variance);
		double base = dataset.getBaselineLOHEstimate();
		double stepSize = 1-base;
		if(variance<=0){
			variance = Math.abs(stepSize);
		}

		
		hmm.setOpdf(0, new OpdfGaussian(1-stepSize,variance ));
		hmm.setOpdf(1, new OpdfGaussian(1,variance ));
		hmm.setOpdf(2, new OpdfGaussian(1+stepSize,variance));
		hmm.setOpdf(3, new OpdfGaussian(1+(2*stepSize),variance));
		
		hmm.setAij(0, 1, 0.1);
		hmm.setAij(0, 0, 0.9);
		hmm.setAij(0, 2, 0);
		hmm.setAij(0, 3, 0);
		
		hmm.setAij(1, 0, 0.05);
		hmm.setAij(1, 1, 0.9);
		hmm.setAij(1, 2, 0.05);
		hmm.setAij(1, 3, 0);
		
		hmm.setAij(2, 0, 0);
		hmm.setAij(2, 1, 0.08);
		hmm.setAij(2, 2, 0.9);
		hmm.setAij(2, 3, 0.2);
		
		hmm.setAij(3, 0, 0);
		hmm.setAij(3, 1, 0);
		hmm.setAij(3, 2, 0.1);
		hmm.setAij(3, 3, 0.9);
		
		return hmm;
	}

	// private static Hmm<ObservationInteger> getHMM_Old(DataSet dataset) {
	//
	// int[] countn = new int[4];
	// List<List<WaveletIF>> plist = dataset.getCapInterval();
	// for(List<WaveletIF> list:plist){
	//
	// for(WaveletIF wi:list){
	//
	// double cn = wi.getCN();
	// int n = (int)(cn/0.5);
	// if(n>=4)n=4;
	// countn[n-1] = countn[n-1]+1;
	//
	// }
	//
	// }
	// int nodesize = 0;
	// int total = 0;
	// for(int m: countn){
	// if(m>0){
	// nodesize++;
	// }
	// total = total+m;
	// }
	// ////
	//
	//
	// OpdfIntegerFactory factory = new OpdfIntegerFactory(4);
	// Hmm<ObservationInteger> hmm = new Hmm<ObservationInteger>(4, factory);
	// int idx = 0;
	// for(int m: countn){
	//
	// double p = (double)m/(double)total;
	// System.out.println(idx+"\t"+p);
	// if(idx>=nodesize)break;
	// hmm.setPi(idx, p);
	// idx++;
	// }
	//
	// hmm.setOpdf(0, new OpdfInteger(new double[] { 0.9, 0.1,0, 0 }));
	// hmm.setOpdf(1, new OpdfInteger(new double[] { 0.05, 0.9,0.05,0 }));
	// hmm.setOpdf(2, new OpdfInteger(new double[] { 0, 0.08,0.9,0.02 }));
	// hmm.setOpdf(3, new OpdfInteger(new double[] { 0, 0,0.1,0.9 }));
	//
	//
	// hmm.setAij(0, 1, 0.1);
	// hmm.setAij(0, 0, 0.9);
	// hmm.setAij(0, 2, 0);
	// hmm.setAij(0, 3, 0);
	//
	// hmm.setAij(1, 0, 0.05);
	// hmm.setAij(1, 1, 0.9);
	// hmm.setAij(1, 2, 0.05);
	// hmm.setAij(1, 3, 0);
	//
	// hmm.setAij(2, 0, 0);
	// hmm.setAij(2, 1, 0.08);
	// hmm.setAij(2, 2, 0.9);
	// hmm.setAij(2, 3, 0.2);
	//
	// hmm.setAij(3, 0, 0);
	// hmm.setAij(3, 1, 0);
	// hmm.setAij(3, 2, 0.1);
	// hmm.setAij(3, 3, 0.9);
	//
	// return hmm;
	// }

}
