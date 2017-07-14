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

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;

public class AdjustBetweenState {

	public static void calc(DataSet dataset) throws IOException {

	
		//
		double baselineestimate = dataset.getBaselineLOHEstimate();
		List<List<WaveletIF>> plist = dataset.getCapInterval();
		for (List<WaveletIF> list : plist) {
			
			double area = 0;
			double totalarea = 0;
			for (WaveletIF wi : list) {
				
				double val = wi.getValue();
				double cn = wi.getHMMValue();
				double diff = getDiff(baselineestimate,cn,val);
				CapInterval ci = (CapInterval)wi;
				area = area+ (diff*ci.getLength());
				totalarea = totalarea +ci.getLength();
			}
			double chromfitting = area/totalarea;
			if(chromfitting>0.5){
				//check between states
				Hmm<ObservationReal>  hmm =
					getHMMWithBetweenStates(list,baselineestimate);
				List<ObservationReal> olist = getList(list);
				int[] hmmary = hmm.mostLikelyStateSequence(olist);
				boolean lessStatesChange = checkStateNo(list,hmmary);
				if(lessStatesChange){
					int cnt = 0;
					for (int m : hmmary) {
						//
						WaveletIF wif = list.get(cnt);
						wif.setHMMValue(states[m]);
						cnt++;

					}
				}
				
			}

		}

	}
	private static boolean checkStateNo(List<WaveletIF> list, int[] hmmary) {
		int stateChangeOrg = getStateChange(list);
		int stateChangeBetween = getStateChange(hmmary);
		
		return stateChangeOrg>=stateChangeBetween;
	}
	
	
	private static int getStateChange(int[] hmmary) {
		int st = 0;
		int cnp =0;
		for(int cn:hmmary){
			
			if(cn!=cnp){
				st++;
			}
			cnp = cn;
		}
		return st;
	}
	private static int getStateChange(List<WaveletIF> list) {
		
		int st = 0;
		double cnp =0;
		for(WaveletIF wi:list){
			double cn = wi.getCN();
			if(cn!=cnp){
				st++;
			}
			cnp = cn;
		}
		return st;
	}


	public static float[] states = new float[]{0.5f,0.75f,1,1.25f,1.5f,1.75f,2};
	
	private static List<ObservationReal> getList(List<WaveletIF> list) {

		List<ObservationReal> olist = new ArrayList<ObservationReal>();
		for (WaveletIF wi : list) {

			ObservationReal oi = new ObservationReal(wi.getDenioseValue());
			olist.add(oi);
		}
		return olist;
	}

	private static double getDiff(double baselineestimate, double cn, double val) {
		
		//
		double scale = 1- baselineestimate;
		if(cn==0.5){
			return Math.abs(baselineestimate-val)/scale;
		}else if(cn==1){
			return Math.abs(1-val)/scale;
		}else if(cn==1.5){
			return Math.abs((2-baselineestimate)-val)/scale;
		}else if(cn==2){
			return Math.abs((3-(2*baselineestimate))-val)/scale;
		}
		return 0;
	}
	
	
	private static Hmm<ObservationReal> getHMMWithBetweenStates(List<WaveletIF> list,double base) {

		int[] countn = new int[10];
		SummaryStatistics ss = new SummaryStatistics();
		Set<Double> checkReg = new HashSet<Double>();
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
		
		int nodesize = 0;
		int total = 0;
		for (int m : countn) {
			if (m > 0) {
				nodesize++;
			}
			total = total + m;
		}
		//asuume normal distribution for each n=1, 1.5, 2, 2.5 , 3, 3.5, 4
		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(7, factory);
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
		double stepSize = (1-base)/2;
		if(variance<=0){
			variance = Math.abs(stepSize);
		}
		System.out.println("variance="+variance);	
		
		
		
		hmm.setOpdf(0, new OpdfGaussian(1-(stepSize*2),variance ));
		hmm.setOpdf(1, new OpdfGaussian(1-(stepSize),variance ));
		hmm.setOpdf(2, new OpdfGaussian(1,variance ));
		hmm.setOpdf(3, new OpdfGaussian(1+stepSize,variance));
		hmm.setOpdf(4, new OpdfGaussian(1+(2*stepSize),variance));
		hmm.setOpdf(5, new OpdfGaussian(1+(3*stepSize),variance ));
		hmm.setOpdf(6, new OpdfGaussian(1+(4*stepSize),variance));

		
		hmm.setAij(0, 1, 0.25);
		hmm.setAij(0, 0, 0.95);
		hmm.setAij(0, 2, 0.25);
		hmm.setAij(0, 3, 0);
		hmm.setAij(0, 4, 0);
		hmm.setAij(0, 5, 0);
		hmm.setAij(0, 6, 0);			
		
		hmm.setAij(1, 0, 0.25);
		hmm.setAij(1, 1, 0.95);
		hmm.setAij(1, 2, 0.25);
		hmm.setAij(1, 3, 0);
		hmm.setAij(1, 4, 0);
		hmm.setAij(1, 5, 0);
		hmm.setAij(1, 6, 0);
				
		hmm.setAij(2, 0, 0.0125);
		hmm.setAij(2, 1, 0.0125);
		hmm.setAij(2, 2, 0.95);
		hmm.setAij(2, 3, 0.0125);
		hmm.setAij(2, 4, 0.0125);
		hmm.setAij(2, 5, 0);
		hmm.setAij(2, 6, 0);
	
		hmm.setAij(3, 0, 0);
		hmm.setAij(3, 1, 0.025);
		hmm.setAij(3, 2, 0.025);
		hmm.setAij(3, 3, 0.95);
		hmm.setAij(3, 4, 0);
		hmm.setAij(3, 5, 0);
		hmm.setAij(3, 6, 0);
		
		hmm.setAij(4, 0, 0);
		hmm.setAij(4, 1, 0);
		hmm.setAij(4, 2, 0.0125);
		hmm.setAij(4, 3, 0.0125);
		hmm.setAij(4, 4, 0.95);
		hmm.setAij(4, 5, 0.0125);
		hmm.setAij(4, 6, 0.0125);
		
		hmm.setAij(5, 0, 0);
		hmm.setAij(5, 1, 0);
		hmm.setAij(5, 2, 0);
		hmm.setAij(5, 3, 0);
		hmm.setAij(5, 4, 0.25);
		hmm.setAij(5, 5, 0.95);
		hmm.setAij(5, 6, 0.25);
		
		hmm.setAij(6, 0, 0);
		hmm.setAij(6, 1, 0);
		hmm.setAij(6, 2, 0);
		hmm.setAij(6, 3, 0);
		hmm.setAij(6, 4, 0.025);
		hmm.setAij(6, 5, 0.025);
		hmm.setAij(6, 6, 0.95);
		
		return hmm;
	}


}
