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

import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.karkinos.ploidy.MatchMatrixBean;

public class PeaksInfo {

	List<Peak> peaklist;
	int[] signalcount;
	double[] ma;
	double[] peaksignals;
	float ploidy =0;
	
	public float getPloidy() {
		return ploidy;
	}

	public void setPloidy(float ploidy) {
		this.ploidy = ploidy;
	}

	List<Double> peakdist;
	List<Peak> mainpeakdist;
	
	MatchMatrixBean matchmatrix;

	public MatchMatrixBean getMatchmatrix() {
		return matchmatrix;
	}

	public void setMatchmatrix(MatchMatrixBean matchmatrix) {
		this.matchmatrix = matchmatrix;
	}

	public List<Peak> getMainpeakdist() {
		return mainpeakdist;
	}

	public void setMainpeakdist(List<Peak> mainpeakdist) {
		this.mainpeakdist = mainpeakdist;
	}

	public double[] getPeaksignals() {
		return peaksignals;
	}

	public void setPeaksignals(double[] peaksignals) {
		this.peaksignals = peaksignals;
	}

	public List<Peak> getPeaklist() {
		return peaklist;
	}

	public void setPeaklist(List<Peak> peaklist) {
		this.peaklist = peaklist;
	}

	public int[] getSignalcount() {
		return signalcount;
	}

	public void setSignalcount(int[] signalcount) {
		this.signalcount = signalcount;
	}

	public double[] getMa() {
		return ma;
	}

	public void setMa(double[] ma) {
		this.ma = ma;
	}

	public double getVal(double x) {

		double sum = 0;
		for (Peak p : peaklist) {

			double r = p.getR();
			double d = getNdistP(x, p.getU(), p.getV());
			sum = sum + (r * d);

		}
		return sum;

	}

	public Number getAcutualVal(double x) {
		double sum = 0;
		for (Peak p : peaklist) {
			if (p.artifitial)
				continue;
			double r = p.getR();
			double d = getNdistP(x, p.getU(), p.getV());
			sum = sum + (r * d);

		}
		return sum;
	}

	public Number getArtifitialVal(double x) {
		double sum = 0;
		for (Peak p : peaklist) {
			if (p.artifitial) {
				double r = p.getR();
				double d = getNdistP(x, p.getU(), p.getV());
				sum = sum + (r * d);
			}

		}
		return sum;
	}

	public static double getNdistP(double x, double u, double v) {

		double p = Math.exp(-0.5 * (Math.pow((x - u), 2) / v))
				/ Math.sqrt(2 * Math.PI * v);
		return p;
	}

	public List<Peak> getCopy() {

		List<Peak> list = new ArrayList<Peak>();
		for (Peak p : peaklist) {
			list.add(p.deepCopy());
		}
		return list;
	}

	public double getSD() {
		double maxsd = 0;
		for (Peak p : peaklist) {
			if(p.getSD()>maxsd){
				maxsd = p.getSD();
			}
		}
		return maxsd;
	}
	

}
