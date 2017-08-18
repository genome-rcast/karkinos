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

public class TheoreticalNodes {

	
	float totalploidy;
	float aAllelePloidy;
	float bAllelePloidy;
	double peakpos;
	float[] imbalance4N = new float[100];
	float[] distTo4N = new float[100];
	
	float[] imbalance2N = new float[100];
	float[] distTo2N = new float[100];
	
	public final static int diploid = 2;
	public final static int tetraploid = 4;
	
	public TheoreticalNodes(int[] ploidy,double peakpos) {
		
		aAllelePloidy = ploidy[0];
		bAllelePloidy = ploidy[1];
		totalploidy = ploidy[0]+ploidy[1];
		this.peakpos = peakpos;
		setUp(diploid,ploidy);
		setUp(tetraploid,ploidy);
		
		
	}
	
	public TheoreticalNodes() {
		// TODO Auto-generated constructor stub
	}

	double diffu;
	double diffpd;
	
	
	public double getDiffu() {
		return diffu;
	}



	public void setDiffu(double diffu) {
		this.diffu = diffu;
	}



	public double getDiffpd() {
		return diffpd;
	}



	public void setDiffpd(double diffpd) {
		this.diffpd = diffpd;
	}



	private void setUp(int i, int[] ploidy) {
		
		
		if(i==diploid){
	    	VirtualPeakCalculator.calc(ploidy,i,imbalance2N,distTo2N,peakpos);	
		}else{
			VirtualPeakCalculator.calc(ploidy,i,imbalance4N,distTo4N,peakpos);				
		}
		
	}



	public double getPeakpos() {
		return peakpos;
	}


	public void setPeakpos(double peakpos) {
		this.peakpos = peakpos;
	}


	public float getTotalploidy() {
		return totalploidy;
	}
	public void setTotalploidy(float totalploidy) {
		this.totalploidy = totalploidy;
	}
	public float getaAllelePloidy() {
		return aAllelePloidy;
	}
	public void setaAllelePloidy(float aAllelePloidy) {
		this.aAllelePloidy = aAllelePloidy;
	}
	public float getbAllelePloidy() {
		return bAllelePloidy;
	}
	public void setbAllelePloidy(float bAllelePloidy) {
		this.bAllelePloidy = bAllelePloidy;
	}
	public float[] getImbalance4N() {
		return imbalance4N;
	}
	public void setImbalance4N(float[] imbalance4n) {
		imbalance4N = imbalance4n;
	}
	public float[] getDistTo4N() {
		return distTo4N;
	}
	public void setDistTo4N(float[] distTo4N) {
		this.distTo4N = distTo4N;
	}
	public float[] getImbalance2N() {
		return imbalance2N;
	}
	public void setImbalance2N(float[] imbalance2n) {
		imbalance2N = imbalance2n;
	}
	public float[] getDistTo2N() {
		return distTo2N;
	}
	public void setDistTo2N(float[] distTo2N) {
		this.distTo2N = distTo2N;
	}


	public String getID() {

		return aAllelePloidy + "-" +bAllelePloidy;
	}
	
}
