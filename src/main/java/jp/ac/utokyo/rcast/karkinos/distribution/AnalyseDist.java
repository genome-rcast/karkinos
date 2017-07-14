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
package jp.ac.utokyo.rcast.karkinos.distribution;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;

public class AnalyseDist implements java.io.Serializable {

	// public SummaryStatistics getTumorratioFromLOH() {
	// return tumorratioFromLOH;
	// }
	//
	// public SummaryStatistics getTumorratioFromGAIN() {
	// return tumorratioFromGAIN;
	// }
	//
	// public SummaryStatistics getTumorratioALL() {
	// return tumorratioALL;
	// }

	public Map<Float, DataHolderByCN> getMap() {
		return map;
	}

	Map<Float, DataHolderByCN> map = new LinkedHashMap<Float, DataHolderByCN>();
	// SummaryStatistics tumorratioFromLOH = new SummaryStatistics();
	// SummaryStatistics tumorratioFromGAIN = new SummaryStatistics();
	// SummaryStatistics tumorratioALL = new SummaryStatistics();
	TumorRatioBean tumorratioFromLOH = null;
	TumorRatioBean tumorratioFromGAIN = null;
	TumorRatioBean tumorratioFromSomatic = null;

	private static final float LOH = 0.5f;
	private static final float twoN = 1f;
	private static final float GAIN = 1.5f;
	private int cntCNV = 0;
	private int lohcount = 0;
	boolean n2havemorepeakdist = false;
	boolean exceedTC = false;
	
	public void analyseDist(DataSet dataset) {

		initmap();
		cntCNV = dataset.getCniVaridateList().size();
		lohcount = 0;
		for(CopyNumberInterval cni:dataset.getCniVaridateList()){
			if(cni.isAllelic())continue;
			if(cni.getCopynumber()==1.0f){
				lohcount++;
			}
		}
		//
		for (SNVHolder snv : dataset.getSnvlist()) {

			float f = (float) snv.getCi().getVaridateVal();
			DataHolderByCN dh = map.get(f);
			//
			if(f==LOH){
				if(overlapWithAllelicGain(snv,dataset.getCniVaridateList())){
					continue;
				}
			}
			if (dh != null) {
				dh.add(snv);
			}
		}
		
		// calculate s.d for tumor n=2;
		DataHolderByCN dh = map.get(1f);
		double mean = dh.tumorhetrosnpMeanSd.getMean();
		double tumorsd = dh.tumorhetrosnpMeanSd.getStandardDeviation();
		tumorratioFromLOH = new TumorRatioBean();
		float[] ret = DistributionFitting.getObserveRatio(map.get(LOH), LOH,
				tumorsd);
		float observedRLOH = ret[0];
		float correl = ret[1];
		if (observedRLOH < 55) {
			// fail to calculata val use val from base line detection
			tumorratioFromLOH.setMode(1);
			observedRLOH = toObserveRatio(dataset.getBaselineLOHEstimate());
			if (observedRLOH < 55) {
				tumorratioFromLOH.setMode(2);
				observedRLOH = 55;
			}

		} else {
			tumorratioFromLOH.setMode(0);
		}
		
		//
		float[] ret2n = DistributionFitting.getObserveRatio(map.get(twoN), LOH,
				tumorsd);
		if(ret[0] < ret2n[0]){
			//2n val have larger peak dist
			n2havemorepeakdist = true;
		}

		tumorratioFromLOH.setObservedratio(observedRLOH);
		int number = getNumber(map.get(LOH).snpDistT);
		tumorratioFromLOH.setNumber(number);
		tumorratioFromLOH.setSd((float) tumorsd);
		float tumorratioLOH = TumorRateCalculator.getTumorRatio(LOH,
				observedRLOH);
		tumorratioFromLOH.setTumorratio(tumorratioLOH);
		tumorratioFromLOH.setCorrel(correl);

		//
		tumorratioFromGAIN = new TumorRatioBean();
		float[] ret0 = DistributionFitting.getObserveRatio(map.get(GAIN), GAIN,
				tumorsd, 50, 67);
		float observedRGAIN = ret0[0];
		float correl0 = ret0[1];
		tumorratioFromGAIN.setCorrel(correl0);
		tumorratioFromGAIN.setObservedratio(observedRGAIN);
		int numberg = getNumber(map.get(GAIN).snpDistT);
		tumorratioFromGAIN.setNumber(numberg);
		tumorratioFromGAIN.setSd((float) tumorsd);
		float tumorratioGAIN = TumorRateCalculator.getTumorRatio(GAIN,
				observedRGAIN);
		tumorratioFromGAIN.setTumorratio(tumorratioGAIN);

	}

	private boolean overlapWithAllelicGain(SNVHolder snv, List<CopyNumberInterval> list) {
		
		//
		if(list==null){
			return false;
		}
		for(CopyNumberInterval cni:list){
			if(cni.isAllelic()&&(cni.getCopynumber()==3)){
				if(cni.getChr().equals(snv.getChr())){
					boolean include 
						= cni.getStart()<=snv.getPos()&& cni.getEnd()>=snv.getPos();
					if(include){
						return true;
					}
				}
			}
		}	
		return false;
	}

	
	public int getTcflg() {
		return tcflg;
	}

	public void setTumorratioFromGAIN(TumorRatioBean tumorratioFromGAIN) {
		this.tumorratioFromGAIN = tumorratioFromGAIN;
	}

	public void reanalyseTC(DataSet dataset) {

		initmap();
		//
		double tcnow = getTumorratio();
		if(tcnow==0){
			tcnow=1;
		}
		int totalpass = 0;
		int exceed1 = 0;
		for (SNVHolder snv : dataset.getSnvlist()) {

			float f = (float) snv.getCi().getVaridateVal();
			DataHolderByCN dh = map.get(f);
			if (dh != null) {
				dh.add(snv);
			}
			if (snv.getFilterResult() != null
					&& snv.getFilterResult().isPassFilter()) {
					totalpass++;
					if((snv.getTumor().getRatio()/tcnow) > 1){
						exceed1++;
					}
			}			

		}
		if(ratiodev(exceed1,totalpass)>0.05){
			exceedTC = true;
		}else{
			exceedTC = false;
		}
		//
		int baseploidy = dataset.getBaseploidy();
		float cn = 1f;
		if(baseploidy == 4){
			cn = 2f;
		}
		DataHolderByCN dh = map.get(cn);
		double mean = dh.tumorhetrosnpMeanSd.getMean();
		double tumorsd = dh.tumorhetrosnpMeanSd.getStandardDeviation();
		tumorratioFromSomatic = new TumorRatioBean();
		float[] rets = DistributionFitting2.getObserveRatio(map.get(cn),
				tumorsd, 5, 50);
		float observedS = rets[0];
		float correls = rets[1];
		tumorratioFromSomatic.setCorrel(correls);
		tumorratioFromSomatic.setObservedratio(observedS);
		int nums = count(map.get(cn).mutationDistTFinalFilter100xdepth);
		tumorratioFromSomatic.setNumber(nums);
		//
		tumorratioFromSomatic.setSd((float) tumorsd);
		float tumorratioSomatic = TumorRateCalculator
				.getTumorRatioSomatic(observedS);
		tumorratioFromSomatic.setTumorratio(tumorratioSomatic);

	}

	// private void addTumorratioStats(SNVHolder snv) {
	//
	// if(!snv.isHetroSNP())return;
	//
	// float f = (float) snv.getCi().getHMMValue();
	// if((f==LOH)||(f==GAIN)){
	//
	// float observedR =snv.getTumor().getRatio();
	// float tumorratio = TumorRateCalculator.getTumorRatio(f, observedR);
	// if(tumorratio==0)return;
	// //20120705-
	// if(underThres(snv.getTumor())||underThres(snv.getNormal()))return;
	// if(f==LOH){
	// tumorratioFromLOH.addValue(tumorratio);
	// }else if(f==GAIN){
	// tumorratioFromGAIN.addValue(tumorratio);
	// }
	// tumorratioALL.addValue(tumorratio);
	// }
	//
	// }

	private double ratiodev(int exceed1, int totalpass) {
		double d = (double)exceed1;
		double d2 = (double)totalpass;
		if(d2==0)return 0;
		return d/d2;
	}

	private int count(float[] dist) {
		// TODO Auto-generated method stub
		try {
			if (dist == null)
				return 0;
			float sum = 0;
			for (float f : dist) {
				sum = sum + f;
			}
			return (int) sum;
		} catch (Exception ex) {

		}
		return 0;
	}

	public TumorRatioBean getTumorratioFromSomatic() {
		return tumorratioFromSomatic;
	}

	private float toObserveRatio(double baselineLOHEstimate) {

		double d = (baselineLOHEstimate - 0.5);
		return (float) ((1 - d) * 100);
	}

	public TumorRatioBean getTumorratioFromLOH() {
		return tumorratioFromLOH;
	}

	public TumorRatioBean getTumorratioFromGAIN() {
		return tumorratioFromGAIN;
	}

	private int getNumber(float[] snpDistT) {
		int sum = 0;
		for (int n = 1; n <= 99; n++) {
			sum = sum + (int) snpDistT[n];
		}
		return sum;
	}

	public static float getGain() {
		return GAIN;
	}

	private boolean underThres(PileUPResult pr) {
		if (pr.getTotalcnt() < (KarkinosProp.mindepth * 2)) {
			return true;
		}
		return false;
	}

	private void initmap() {

		for (float f = 0.5f; f <= 2; f = f + 0.25f) {
			map.put(f, new DataHolderByCN());
		}

	}


	int tcflg = 0;
	public float getTumorratio() {

		float tr1 = tumorratioFromLOH.getTumorratio();
		float tr2 = tumorratioFromGAIN.getTumorratio();

		float tr0 = 0f;

		if (tumorratioFromSomatic != null) {
			tr0 = tumorratioFromSomatic.getTumorratio();
		}
		float ret = 0f;
		boolean toomuchcnv = lohcount >= 50;
		if (!toomuchcnv && (tr1 > 0) && (tumorratioFromLOH.getNumber() >= 50)) {
			ret = tr1;
			tcflg = 1;

		} else if (!toomuchcnv && (tr2 > 0) && (tumorratioFromGAIN.getNumber() >= 50)) {
			ret = tr2;
			tcflg = 2;

		} else {
			ret = tr0;
			tcflg = 3;

		}
		
		if (tcflg == 2) {
			if (tr0 > 0 && tr0 < 0.5 && (ret > tr0) && (ret - tr0) > 0.3) {
				ret = Math.min(tr1, tr2);
				if ((ret > tr0) && (ret - tr0) > 0.3) {
					ret = tr0;
					tcflg = 3;
				}
			}
		}
		if (tcflg == 1) {
			if (tumorratioFromSomatic != null && n2havemorepeakdist) {
				tr0 = tumorratioFromSomatic.getTumorratio();
				if(tumorratioFromSomatic.number>=40){
					
					if(tr1+0.25 < tr0){
					 ret = tr0;
					 tcflg = 3;
					}
					
				}
			}

		}		
		if(exceedTC){
			if(ret<tr0){
			  ret = tr0;
			  tcflg = 3;
			} 
		}
		return ret;

		// if(tumorratioFromLOH.getNumber()>tumorratioFromGAIN.getCorrel()){
		// return tr1;
		// }
		// if(tumorratioFromLOH.getCorrel()>tumorratioFromGAIN.getCorrel()){
		// return tr1;
		// }else{
		// return tr2;
		// }

	}

}
