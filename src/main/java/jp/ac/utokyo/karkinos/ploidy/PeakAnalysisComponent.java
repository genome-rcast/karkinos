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
package jp.ac.utokyo.karkinos.ploidy;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.SNVHolderPlusACnv;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;

public class PeakAnalysisComponent {

	public float getHpeakdistance() {
		return hpeakdistance;
	}



	AFCounter normalSNP = new AFCounter();
	AFCounter tumorSNP = new AFCounter();
	AFCounter somaticSNV = new AFCounter();
	Map<String,ChrAllelicPeak> ap = new LinkedHashMap<String,ChrAllelicPeak>();
	
	
	AFCounter tumorSNPEven = new AFCounter();
	AFCounter tumorSNPOdd = new AFCounter();
	
	boolean evenpeak = true;
	float hpeakdistance = 0;
	float aAfreq = 0;
	float bAfreq = 0;
	int numdistpeak =0;
	boolean complexpeak = false;
	double sdnormal = 0;

	public void setSNPInfo(SNVHolder snv) {

		double normalAF = snv.getNormal().getRatio();
		normalSNP.setValue(normalAF);
		double tumorAF = snv.getTumor().getRatio();
		tumorSNP.setValue(tumorAF);

	}

	public void setSMutaion(SNVHolder snv) {

		double tumorAF = snv.getTumor().getRatio();
		somaticSNV.setValue(tumorAF);

	}

	public void setSNPPlus(SNVHolderPlusACnv snvHolderPlusACnv) {

		if (snvHolderPlusACnv != null) {
			
			//
			String chr = snvHolderPlusACnv.getSnv().getChr();
			ChrAllelicPeak cap = null;
			if(ap.containsKey(chr)){
				cap = ap.get(chr);
			}else{
				cap = new ChrAllelicPeak();
				ap.put(chr, cap);
			}
			cap.add(snvHolderPlusACnv);
			
		}

	}
	
	Set<String> evenChrSet = new HashSet<String>();
	Set<String> oddChrSet = new HashSet<String>();
	
	private static final float eventhres = 0.525f;
	public void analyse() {

		//
		sdnormal = normalSNP.getSs().getStandardDeviation();
		Set<String> set = ap.keySet();
		Iterator<String> ite = set.iterator();
		//
		Set<Boolean> sb = new HashSet<Boolean>();
		int evencnt = 0;
		int oddcnt = 0;
		
		while(ite.hasNext()){
			
			String chr = ite.next();
			if(chr.toUpperCase().contains("X")||chr.toUpperCase().contains("Y")){
				continue;
			}			
			ChrAllelicPeak cap = ap.get(chr);
			float pd = cap.getTumorSNP().getPeakDistance(sdnormal)[0];			
			int np = (int)cap.getTumorSNP().getPeakDistance(sdnormal)[1];
			System.out.println(pd+"\t"+chr);
			boolean even = (pd < eventhres) && (np <= 1);
			//
			cap.setEven(even);
			sb.add(even);		
			if(even){
				evenChrSet.add(chr);
			}else{
				oddChrSet.add(chr);
			}
			

			for(SNVHolderPlusACnv snv:cap.list){
				double d = snv.getSnv().getTumor().getRatio();
				if(even){
					
					//					
					if(middle(d)){
						tumorSNPEven.setValue(d);
						evencnt++;
					}else{
						tumorSNPOdd.setValue(d);
						oddcnt++;
					}
					
				}else{
					
					//
					tumorSNPOdd.setValue(d);
					oddcnt++;
				}
			}
			
		
			
		}		
		if(sb.size()>1){
			complexpeak = true;
			evenpd = tumorSNPEven.getPeakDistance(sdnormal)[0];
			oddpd = tumorSNPOdd.getPeakDistance(sdnormal)[0];
			evenr = ((double)evencnt/(double)(evencnt+oddcnt));
			oddr = ((double)oddcnt/(double)(evencnt+oddcnt));
		}
		
		//
		hpeakdistance = tumorSNP.getPeakDistance(sdnormal)[0];
		numdistpeak = (int)tumorSNP.getPeakDistance(sdnormal)[1];
		evenpeak = (hpeakdistance < eventhres) && (numdistpeak <= 1);		
		
		
		
	}
	///
	
	
	
	private boolean middle(double d) {
		if(d<0.55 && d > 0.45){
			return true;
		}
		return false;
	}

	public Set<String> getEvenChrSet() {
		return evenChrSet;
	}

	public Set<String> getOddChrSet() {
		return oddChrSet;
	}

	double evenpd;
	public double getEvenPD() {
		// 
		return evenpd;
	}
	
	double evenr;
	public double getEvenR() {
		// 
		return evenr;
	}
	
	double oddr;
	public double getOddR() {
		// 
		return oddr;
	}
	
	double oddpd;
	public double getOddPD() {
		// 
		return oddpd;
	}

	public boolean isComplexPeak() {
		
		return complexpeak;
	}



}
