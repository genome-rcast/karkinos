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
package jp.ac.utokyo.rcast.karkinos.summary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class ProbabilityUtils {

	public static void main(String[] arg) {

		int cnt0 = 0;
		int cnt1 = 0;

		int n = 50000;

		double p = 3.2391456991570703E-6;
		int len = 500;

		Date start = Calendar.getInstance().getTime();
//		for (int m = 0; m < n; m++) {
//
//			int s = getBinomialP(len, p);
//			if (s == 0) {
//				cnt0++;
//			} else if (s == 1) {
//				cnt1++;
//			}
//
//		}
		System.out.println(getPval(p, 1,len));
		
		Date end = Calendar.getInstance().getTime();
		System.out.println("cnt0=" + cnt0 + " cnt1=" + cnt1);
		System.out.println("time=" + (end.getTime() - start.getTime()));
	}

	static BinomialDistribution bdist = new BinomialDistributionImpl(0, 0);
	public static double getPval(double p, int cnt, int genelen) {
		
		bdist.setNumberOfTrials(genelen);
		bdist.setProbabilityOfSuccess(p);
		//BinomialDistribution bdist = new BinomialDistributionImpl(genelen, p);
		return bdist.probability(cnt);

	}

	static final int similationtrialn = 10000;
	static RandomDataImpl rd = new RandomDataImpl();

	public static double getPValbySample(
			Map<String, Double> backGroundMutationRate, int genelength,
			double likehood) throws MathException {

		Map<String, BinomialDistribution> distmap = new HashMap<String, BinomialDistribution>();

		for (String sample : backGroundMutationRate.keySet()) {

			double p = backGroundMutationRate.get(sample);
			BinomialDistribution bdist = new BinomialDistributionImpl(
					genelength, p);
			distmap.put(sample, bdist);

		}
		//
		TreeMap<Double,MyCounter> simdata 
			= new TreeMap<Double,MyCounter>();
		
		int simsize = similationtrialn;

		for (int n = 0; n < similationtrialn; n++) {
			double testlikehood = getlh(distmap, genelength);
			if(simdata.containsKey(testlikehood)){
				simdata.get(testlikehood).inc();
			}else{
				simdata.put(testlikehood, new MyCounter());
			}
					
		}
		double firstval = simdata.firstEntry().getKey();

		if (likehood > firstval) {
			
			simsize = simsize + (similationtrialn*9);
			for (int n = 0; n < (similationtrialn*9); n++) {
				double testlikehood = getlh(distmap, genelength);
				if(simdata.containsKey(testlikehood)){
					simdata.get(testlikehood).inc();
				}else{
					simdata.put(testlikehood, new MyCounter());
				}
				
			}
		}	
		
		//
							
		int idx = 0;
		int mid = 0;
		if (simdata.containsKey(likehood)) {
			
			int cnt = simdata.get(likehood).getCnt();
			mid = cnt/2;			
		}
		idx = getIndex(simdata,likehood)+mid;
		//				
		double pval =  (double) (idx + 1) / (double) (simsize + 1);
		return pval;

	}

	private static int getIndex(TreeMap<Double, MyCounter> simdata,
			double likehood) {
		
		//
		SortedMap<Double, MyCounter> map  = simdata.tailMap(likehood);
		if(map==null){
			return 0;
		}
		
		Iterator<Double> ite = map.keySet().iterator();
		int idx = 0;
		while(ite.hasNext()){
			
			idx = idx + map.get(ite.next()).getCnt();
			
		}
		return idx;
	}

	private static double getlh(Map<String, BinomialDistribution> distmap,
			int genelen) {

		double likehood = 0;
		//
		Iterator<String> keys = distmap.keySet().iterator();

		while (keys.hasNext()) {

			String key = keys.next();
			BinomialDistribution bdist = distmap.get(key);
			int randcnt = getBinomialP(genelen, bdist.getProbabilityOfSuccess());
			double pvaleach = bdist.probability(randcnt);
			likehood = likehood + (-1 * Math.log(pvaleach));
		}
		return likehood;
	}

	public static int getBinomial(int n, double p) {

		try {
			return rd.nextBinomial(n, p);
		} catch (MathException e) {

			return getBinomialOld(n, p);
		}

	}

	public static int getBinomialP(int n, double p) {

		double mu = n * p;
		return (int) rd.nextPoisson(mu);

	}

	public static int getBinomialOld(int n, double p) {

		int x = 0;
		for (int i = 0; i < n; i++) {
			if (Math.random() < p)
				x++;
		}
		return x;

	}

}
