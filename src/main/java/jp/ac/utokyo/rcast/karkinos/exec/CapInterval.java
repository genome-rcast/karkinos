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
package jp.ac.utokyo.rcast.karkinos.exec;


import htsjdk.samtools.SAMRecord;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class CapInterval implements WaveletIF, java.io.Serializable {
	private static float divide(final int i, final int i2) {
		return (float)((double)i/(double)i2);
	}

	private final String chr;
	private int start, end; // 1-based closed range
	private final boolean gene;
	private final int length; // NB: This value is invalid if CapInterval is merged.
	private int totallength;
	private float cgParcent;
	private float duality;
	private int peakIdx;
	private float cnvtotal;
	private float aafreq;
	private float bafreq;
	private float normalAveDepth;
	private double hmmvalue;
	private double varidateval;
	private CNVInfo cnvinfo;
	private boolean startChrom;
	private boolean endChrom;

	public CapInterval(final String chr, final int start, final int end, final boolean gene) {
		this(chr, start, end, gene, 0f, 1f);
	}

	public CapInterval(final String chr, final int start, final int end, final boolean gene,
			final float cgParcent, final float duality) {
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.gene = gene;
		this.length = end - start + 1;
		this.totallength = this.length;
		this.cgParcent = cgParcent;
		this.duality = duality;
		this.peakIdx = 0;
		this.cnvtotal = 0f;
		this.aafreq = 0f;
		this.bafreq = 0f;
		this.normalAveDepth = 0f;
		this.hmmvalue = 0;
		this.varidateval = 0;
		this.cnvinfo = null;
		this.startChrom = false;
		this.endChrom = false;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public boolean isGene() {
		return gene;
	}

	public float getNormalAveDepth() {
		return normalAveDepth;
	}

	public void setNormalAveDepth(final float normalAveDepth) {
		this.normalAveDepth = normalAveDepth;
	}

	public int getPeakIdx() {
		return peakIdx;
	}

	public void setPeakIdx(final int peakIdx) {
		this.peakIdx = peakIdx;
	}

	public float getCnvtotal() {
		return cnvtotal;
	}

	public void setCnvtotal(final float cnvtotal) {
		this.cnvtotal = cnvtotal;
	}

	public float getAafreq() {
		return aafreq;
	}

	public void setAafreq(final float aafreq) {
		this.aafreq = aafreq;
	}

	public float getBafreq() {
		return bafreq;
	}

	public void setBafreq(final float bafreq) {
		this.bafreq = bafreq;
	}

	public boolean isStartChrom() {
		return startChrom;
	}

	public void setStartChrom(final boolean startChrom) {
		this.startChrom = startChrom;
	}

	public boolean isEndChrom() {
		return endChrom;
	}

	public void setEndChrom(final boolean endChrom) {
		this.endChrom = endChrom;
	}

	public int getLength() {
		return length;
	}

	public float getCgParcent() {
		return cgParcent;
	}

	public String getInfoStr() {
		StringBuffer sb = new StringBuffer();
		sb.append(chr+"\t");
		sb.append(start+"\t");
		sb.append(end+"\t");
		sb.append(cnvinfo.getNormalcnt()+"\t");
		sb.append(cnvinfo.getTumorcnt()+"\t");
		sb.append(cnvinfo.getTnratio()+"\t");
		sb.append(cnvinfo.getDenoise()+"\t");
		sb.append(cnvinfo.getCopynumber()+"\t");
		sb.append(getHMMValue()+"\t");
		sb.append(getVaridateVal());
		return sb.toString();
	}

	public void setCgParcent(final float cgParcent) {
		this.cgParcent = cgParcent;
	}

	public float getDepth(final long totalbase) {
		final double len = end - start + 1;
		final double dp = (double) totalbase / len;
		return (float) dp;
	}

	public void merge(final CapInterval iv) {
		start = Math.min(start, iv.start);
		end = Math.max(end, iv.end);
		totallength = totallength + iv.length;
		duality = divide(totallength, end - start + 1);
//		if(cnvinfo==null){
//			cnvinfo = iv.cnvinfo;
//		}else{
//			cnvinfo.merge(iv.cnvinfo);
//		}
	}

	public boolean intersect(final CapInterval iv) {
		return intersect(iv.start, iv.end);
	}

	public boolean intersect(final SAMRecord sr) {
		final int s = sr.getAlignmentStart();
		int e = sr.getAlignmentEnd();
		if (e == 0 || e == s) {
			e = sr.getAlignmentStart() + sr.getReadLength() - 1;
		}
		return intersect(s, e);
	}

	// The arguments take closed range [s,e]
	public boolean intersect(final int s, final int e) {
		return start <= e && s <= end;
	}

	public boolean equals(Object obj) {
		if (!(obj instanceof CapInterval)) {
			return false;
		}
		CapInterval ci0 = (CapInterval)obj;
		return ci0.chr.equals(chr)
			&& (ci0.start == start)
			 && (ci0.end == end);
	}

	public int hashCode() {
		return this.infoStr().hashCode();
	}

	public String infoStr() {
		return chr+":"+start+"-"+end+" ("+(length)+")";
	}

	public void setCNVInfo(final CNVInfo cnvinfo) {
		this.cnvinfo = cnvinfo;
	}

	public CNVInfo getCNVInfo() {
		return cnvinfo;
	}

	public void setDenioseValue(final double denoise) {
		cnvinfo.setDenioseValue(denoise);
	}

	public void setCN(final double copynumber) {
		cnvinfo.setCN(copynumber);
	}

	public double getOriginalValue() {
		return cnvinfo.getOriginalTnratio();
	}

	public double getValue() {
		return cnvinfo.getTnratio();
	}

	public boolean include(final int pos) {
		return include(pos, 0);
	}

	public boolean include(final int pos, final int mergin) {
		return start - mergin <= pos && pos <= end + mergin;
	}

	public double getDenioseValue() {
		return cnvinfo.getDenoise();
	}

	public double getCN() {
		return cnvinfo.getCopynumber();
	}

	public void setHMMValue(final double d) {
		 hmmvalue = d;
		 varidateval = d;
	}

	public double getHMMValue() {
		return hmmvalue;
	}

	public float getDuality() {
		return duality;
	}

	public void setGCAdjustedTNratio(final double adjustedY) {
		if (cnvinfo != null) {
			cnvinfo.tnratio = adjustedY;
		}
	}

	public double getVaridateVal() {
		return varidateval;
	}

	public void setVaridateval(final double varidateval) {
		this.varidateval = varidateval;
	}

}
