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
import java.util.Collections;
import java.util.List;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class BinData {

	int udepth;
	int ldepth;
	SummaryStatistics hetroSNPAF = new SummaryStatistics();
	SummaryStatistics secondAF = new SummaryStatistics();
	//
	List<Point2D> afdepth = new ArrayList<Point2D>();

	public void setHetroSNPAF(int depth, double tr) {
		hetroSNPAF.addValue(tr);
	}

	public void setSecondAF(int depth, double secondAlleleF) {
		secondAF.addValue(secondAlleleF);		
	}
	
	public void setSNV(int depth, double tr) {
		afdepth.add(new Point2D.Double(tr, depth));
	}

	public BinData(int ldepth, int udepth) {

		this.ldepth = ldepth;
		this.udepth = udepth;

	}

	public int getDepth() {
		return (ldepth + udepth) / 2;
	}

	public int getUdepth() {
		return udepth;
	}

	public void setUdepth(int udepth) {
		this.udepth = udepth;
	}

	public int getLdepth() {
		return ldepth;
	}

	public void setLdepth(int ldepth) {
		this.ldepth = ldepth;
	}

	public SummaryStatistics getHetroSNPAF() {
		return hetroSNPAF;
	}

	public void setHetroSNPAF(SummaryStatistics hetroSNPAF) {
		this.hetroSNPAF = hetroSNPAF;
	}

	int numCandidate;
	boolean emExcuted = false;
	double borerAF = 0.15f;

	public void em() {


		try {
			EMCGM emcgm = new EMCGM();
			emcgm.setInitVariance(hetroSNPAF.getVariance());
			emcgm.analyse(afdepth);
			numCandidate = emcgm.getNumcandidate();
			borerAF = emcgm.getBorderAF();
				
			if (emcgm.isSuccess()) {
				emExcuted = true;
				usedefalt = false;
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public List<Point2D> getAfdepth() {
		return afdepth;
	}

	public double getAFBorder() {
		
		double sd = hetroSNPAF.getStandardDeviation();
		double defret = sd*2;
		if (Double.isNaN(sd)) {
			sd = 0;
		}
		if (usedefalt) {
			// case all EM fail, return s.d. of hetro SNP AF diff

			// 3sd
			return defret;
		} else if (!emExcuted) {
			
			if(getSamplesize()<=10){
				return defret;
			}
			// case EM fail for this bin but succeed in at least one other bin
			double r = getUpperAFNoise(predictedCandnum);
			if((r<0)||(r<(defret))){
				return defret;
			}
			return r;
		} else {
			
			System.out.println("EM used="+getDepth()+"\t"+borerAF);
			if((borerAF<0.15) && (getDepth()<40)){
				return defret;
			}			
			return borerAF;
			
		}

	}

	public int getNumCandidate() {
		return numCandidate;
	}

	public double getUpperAFNoise(int candnum) {

		int idx = afdepth.size() - (candnum+1);
		if (idx < 0){
			return 0;
			//no noise
		}	
		if(idx>=afdepth.size()){
			return 0;
		}
		return afdepth.get(idx).getX();

	}

	public void sortList() {

		sort(afdepth);

	}

	private void sort(List<Point2D> afdepth) {

		Collections.sort(afdepth, new MyCompEM());
	}

	public boolean isEmExcuted() {
		return emExcuted;
	}

	boolean usedefalt = true;

	int predictedCandnum = 0;

	public void setPredictedCandnum(int predictedCandnum) {

		this.predictedCandnum = predictedCandnum;
		usedefalt = false;
	}

	public boolean includepeth(int d) {
		
		if(udepth==0){
			udepth = 10000;
		}
		return d>= ldepth && d<=udepth;		
	}

	public int getSamplesize() {
		if(afdepth==null)return 0;
		return afdepth.size();
	}



}
