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

import jp.ac.utokyo.karkinos.noisefilter.NDist;
import jp.ac.utokyo.karkinos.ploidy.PeakAnalysisComponent;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;

public class Peak {

	double x;
	double y;

	double u;
	double v;
	double r;
	double vinit;
	boolean artifitial;
	int peakidx;
	float copynum;

	public float getCopynum() {
		return copynum;
	}

	PeakAnalysisComponent pac;

	public PeakAnalysisComponent getPac() {
		return pac;
	}

	public void setPac(PeakAnalysisComponent pac) {
		this.pac = pac;
	}

	public int getPeakidx() {
		return peakidx;
	}

	public void setPeakidx(int peakidx) {
		this.peakidx = peakidx;
	}

	public boolean isArtifitial() {
		return artifitial;
	}

	public void setArtifitial(boolean artifitial) {
		this.artifitial = artifitial;
	}

	public double getVinit() {
		return vinit;
	}

	public void setVinit(double vinit) {
		this.vinit = vinit;
	}

	double SUMZ;
	double SUMZX;
	double SUMZX2;

	float cn;

	public float getCn() {
		return cn;
	}

	public void setCn(float cn) {
		this.cn = cn;
	}

	public double getSUMZ() {
		return SUMZ;
	}

	public void setSUMZ(double sUMZ) {
		SUMZ = sUMZ;
	}

	public double getSUMZX() {
		return SUMZX;
	}

	public void setSUMZX(double sUMZX) {
		SUMZX = sUMZX;
	}

	public double getSUMZX2() {
		return SUMZX2;
	}

	public void setSUMZX2(double sUMZX2) {
		SUMZX2 = sUMZX2;
	}

	public Peak(double x2, double y2) {

		x = x2;
		y = y2;
	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double getU() {
		return u;
	}

	public void setU(double u) {
		this.u = u;
	}

	public double getV() {
		return v;
	}

	public void setV(double v) {
		this.v = v;
	}

	public double getR() {
		return r;
	}

	public void setR(double r) {
		this.r = r;
	}

	public double getSD() {

		if (v <= 0)
			return 0;
		return Math.sqrt(v);
	}

	boolean defualt;

	public void setDefault(boolean b) {
		defualt = b;
	}

	public boolean getDefault() {
		return defualt;
	}

	public Peak deepCopy() {

		Peak p = new Peak(x, y);

		p.u = u;
		p.v = v;
		p.r = r;

		p.cn = cn;
		p.defualt = defualt;
		p.artifitial = artifitial;
		p.peakidx = peakidx;
		return p;

	}

	public void setXYifYBigger(double x2, double y2) {

		if (y2 > y) {
			y = y2;
			x = x2;
			u = x2;
		}

	}

	List<ChildPeak> childPeakList = null;

	public float setCN(CapInterval ci) {

		float aafreq = 0f;
		float bafreq = 0f;
		float cnvtotal = 0f;
		cnvtotal = 2f;
		ci.setAafreq(1f);
		ci.setBafreq(1f);
		ci.setCnvtotal(cnvtotal);
		ci.setHMMValue((double) cnvtotal / (double) 2);
		int idx = 0;
		if (pac == null) {


			copynum = cnvtotal;
			
			try {
				ChildPeak cp = getChildPeaks().get(idx);
				aafreq = cp.getAaf();
				bafreq = cp.getBaf();
				cnvtotal = cp.getTaf();
				ci.setAafreq(aafreq);
				ci.setBafreq(bafreq);
				ci.setCnvtotal(cnvtotal);
				ci.setHMMValue((double) cnvtotal / (double) 2);
				
				
			} catch (Exception ex) {

			}
			return cnvtotal;

		}else if (pac.isComplexPeak()) {

			int maxidx = 0;
			double maxr = 0;
			for (ChildPeak cp : childPeakList) {

				// ci.getChr();
				if (cp.getChrom() != null) {
					if (cp.getChrom().contains(ci.getChr())) {
						break;
					}
				}
				if (cp.r > maxr) {
					maxr = cp.r;
					maxidx = idx;
				}
				idx++;
			}
			if (idx == getChildPeaks().size()) {
				idx = maxidx;
			}

		}

		ChildPeak cp = getChildPeaks().get(idx);
		aafreq = cp.getAaf();
		bafreq = cp.getBaf();
		cnvtotal = cp.getTaf();
		ci.setAafreq(aafreq);
		ci.setBafreq(bafreq);
		ci.setCnvtotal(cnvtotal);
		ci.setHMMValue((double) cnvtotal / (double) 2);
		copynum = cnvtotal;
		return cnvtotal;
	}

	public List<ChildPeak> getChildPeaks() {

		if (childPeakList == null) {
			childPeakList = new ArrayList<ChildPeak>();
			if (pac == null)
				return childPeakList;

			if (pac.isComplexPeak()) {

				ChildPeak even = new ChildPeak();
				ChildPeak odd = new ChildPeak();
				even.r = pac.getEvenR() * r;
				even.u = u;
				even.peakdist = pac.getEvenPD();
				even.chrom = pac.getEvenChrSet();

				odd.r = pac.getOddR() * r;
				odd.u = u;
				odd.peakdist = pac.getOddPD();
				odd.chrom = pac.getOddChrSet();

				childPeakList.add(even);
				childPeakList.add(odd);

			} else {

				ChildPeak peak = new ChildPeak();
				peak.r = r;
				peak.u = u;
				peak.peakdist = pac.getHpeakdistance();
				childPeakList.add(peak);
			}

		}

		return childPeakList;
	}

}
