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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class HMMCNVAnalysisFromEM {
	
	public static void calc(DataSet dataset, PeaksInfo pi) {

			
		List<List<WaveletIF>> plist = dataset.getCapInterval();
		List<Peak> peaklist = pi.getCopy();
		int chrinx = 1;
		for (List<WaveletIF> list : plist) {

			
			List<Peak> usepeaklist = checkExsistance(peaklist,list);
			Map<Integer,Integer> plistidx = toPeakListIdx(usepeaklist);
			System.out.println(chrinx+"\t"+plistidx);
	
			chrinx++;
			if((usepeaklist.size()==0)||usepeaklist.size()>15){
				setdef(list);
				continue;
			}					
			// create HMM for each chromosome
			Hmm<ObservationReal> hmm = getHMM(usepeaklist);
			//
			List<ObservationReal> olist = getList(list,usepeaklist);
			int[] hmmary = hmm.mostLikelyStateSequence(olist);
						
			int cnt = 0;
			int mb4 = -1;
			for (int m : hmmary) {

				//
				
				mb4 = m;
				WaveletIF wif = list.get(cnt);
				CapInterval ci = (CapInterval)wif;
				ci.setPeakIdx(plistidx.get(m));
				cnt++;
				
				System.out.println(cnt + "\t" + wif.getDenioseValue() + "\t" + m +"\t" + plistidx.get(m) );

			}
			
		}

	}

	private static List<ObservationReal> trim(List<ObservationReal> olist, int i) {
		
		List<ObservationReal> l = new ArrayList<ObservationReal>();
		for(int n = 0;n<i;n++){
			l.add(olist.get(n));
		}
		return l;
	}

	private static void setdef(List<WaveletIF> list) {
		
		for(WaveletIF wif: list){
			wif.setHMMValue(1.0f);
		}
		
	}

	private static Map<Integer, Integer> toPeakListIdx(List<Peak> usepeaklist) {
		
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		int idx = 0;
		for(Peak p:usepeaklist){
			
			map.put(idx, p.getPeakidx());			
			idx++;
			
		}
		return map;
		
	}

	private static List<Peak> checkExsistance(List<Peak> peaklist,
			List<WaveletIF> list) {
		
		List<Peak> peaklistn = new ArrayList<Peak>();
		for(Peak p:peaklist){
			
			boolean c = contain(p,list);
			if(c){
				peaklistn.add(p);
			}
						
		}		
		if(peaklistn.size()==0){
			return peaklist;
		}
		return peaklistn;
	}

	private static boolean contain(Peak p, List<WaveletIF> list) {
		
		
		int total = list.size();
		if(total==0)return false;
		int cin = 0;
		for(WaveletIF wif:list){
			double val = wif.getDenioseValue();
			double sd = p.getSD();
			double s = p.getU()- (sd*3);
			double e = p.getU() + (sd*3);
			if((val>s)&&(val<e)){
				cin++;
			}
		
		}
		double r = (double)((double)cin/(double)total);
		return r > 0.001;
	}
	


	private static List<ObservationReal> getList(List<WaveletIF> list, List<Peak> usepeaklist) {

		
		double min=10,max = -10;
		for(Peak p:usepeaklist){
			double u = p.getU();
			if(min > u){
				min = u;
			}
			if(max < u){
				max = u;
			}			
			
		}
		
		List<ObservationReal> olist = new ArrayList<ObservationReal>();
		for (WaveletIF wi : list) {
			
			double d = wi.getDenioseValue();
			if(d>max+0.5){
				d = max;
			}
			if(d<min-0.5){
				d=min;
			}
			ObservationReal oi = new ObservationReal(d);
			olist.add(oi);
		}
		return olist;
	}

	private static Hmm<ObservationReal> getHMM(List<Peak> peaks) {

		// asuume normal distribution for each n=1,2,3,4
		OpdfGaussianFactory factory = new OpdfGaussianFactory();
		Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(peaks.size(),
				factory);
		//
		int n = 0;
		int idxn0 = 0;
		int idxn1 = 0;
		int idxn2 = 0;
		int idxn3 = 0;
		for (Peak peak : peaks) {

			hmm.setOpdf(n, new OpdfGaussian(peak.getU(), peak.getV()));
			if(peak.getCn()==0.5f){
				idxn0 = n;
			}
			if(peak.getCn()==1.0f){
				idxn1 = n;
			}
			if(peak.getCn()==1.5f){
				idxn2 = n;
			}
			if(peak.getCn()==2.0f){
				idxn3 = n;
			}			
			n++;
		}
		
		double portion = 0.1/(double)(peaks.size()-1);
		for(int l= 0;l<peaks.size();l++){
			
			for(int m= 0;m<peaks.size();m++){
				
				//2n if
				if(l==m){
					
					hmm.setAij(l, m, 0.9);
				
				}else{
														
					hmm.setAij(l, m, portion);
				}
				
			}
			
		}

		return hmm;
	}


	

}
