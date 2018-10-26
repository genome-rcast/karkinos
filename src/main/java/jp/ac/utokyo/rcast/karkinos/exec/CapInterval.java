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

public class CapInterval implements WaveletIF, java.io.Serializable{

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
		return this.geneName != null && !this.geneName.equals("rna");
	}

	public boolean isHLA() {
		return this.geneName != null && this.geneName.contains("HLA");
	}

	public String getGeneName() {
		return this.geneName != null ? this.geneName : "rna";
	}

	String chr;
	int start = 0;
	int end = 0;
	private String geneName;
	long total = 0;
	public float getNormalAveDepth() {
		return normalAveDepth;
	}

	public void setNormalAveDepth(float normalAveDepth) {
		this.normalAveDepth = normalAveDepth;
	}

	int length;
	float cgParcent = 0f;
	float duality = 1f;
	boolean startChrom;
	boolean endChrom;
	
	int peakIdx=0;
	float cnvtotal=0f;
	float aafreq=0f;
	float bafreq=0f;
	float normalAveDepth=0f;

	public int getPeakIdx() {
		return peakIdx;
	}

	public void setPeakIdx(int peakIdx) {
		this.peakIdx = peakIdx;
	}

	public float getCnvtotal() {
		return cnvtotal;
	}

	public void setCnvtotal(float cnvtotal) {
		this.cnvtotal = cnvtotal;
	}

	public float getAafreq() {
		return aafreq;
	}

	public void setAafreq(float aafreq) {
		this.aafreq = aafreq;
	}

	public float getBafreq() {
		return bafreq;
	}

	public void setBafreq(float bafreq) {
		this.bafreq = bafreq;
	}

	public boolean isStartChrom() {
		return startChrom;
	}

	public void setStartChrom(boolean startChrom) {
		this.startChrom = startChrom;
	}

	public boolean isEndChrom() {
		return endChrom;
	}

	public void setEndChrom(boolean endChrom) {
		this.endChrom = endChrom;
	}

	public int getLength() {
		return length;
	}

	public float getCgParcent() {
		return cgParcent;
	}

	int totallength =0;

	public String getInfoStr() {
		
		StringBuffer sb = new StringBuffer();
		sb.append(chr+"\t");
		sb.append(start+"\t");
		sb.append(end+"\t");
		sb.append(cnvinfo.normalcnt+"\t");
		sb.append(cnvinfo.tumorcnt+"\t");
		sb.append(cnvinfo.tnratio+"\t");
		sb.append(cnvinfo.getDenoise()+"\t");
		sb.append(cnvinfo.getCopynumber()+"\t");
		sb.append(getHMMValue()+"\t");
		sb.append(getVaridateVal());
		return sb.toString();
	}
	
	public CapInterval(String _chr, int _start, int _end, String _geneName, float cgp, float _duality) {
		this(_chr,_start,_end,_geneName);
		cgParcent = cgp;
		duality = _duality;
	}

	public void setCgParcent(float cgParcent) {
		this.cgParcent = cgParcent;
	}

	public float getDepth(long totalbase) {

		double len = end - start;
		double dp = (double) totalbase / len;
		return (float) dp;
	}


	public CapInterval(String _chr, int _start, int _end, String _geneName) {
		chr = _chr;
		start = _start;
		end = _end;
		geneName = _geneName;
		length = end-start;
		totallength = length;
	}

	public void merge(CapInterval iv) {
		int s = iv.start;
		int e = iv.end;
		if (s < start) {
			start = s;
		}
		if (end < e) {
			end = e;
		}
		if (!this.isGene() && iv.isGene()) {
			this.geneName = iv.geneName;
		}
		totallength = totallength +iv.length;
		duality = devide(totallength,end-start);
//		if(cnvinfo==null){
//			cnvinfo = iv.cnvinfo;
//		}else{
//			cnvinfo.merge(iv.cnvinfo);
//		}
//		
	}

	private float devide(int i, int i2) {
		
		return (float)((double)(double)i/(double)i2);
	}

	public boolean intersect(CapInterval iv) {
		int s = iv.start;
		int e = iv.end;
		boolean overlap = start <= e && s <= end;
		return overlap;
	}
	
	public boolean intersect(int s,int e) {
		boolean overlap = start <= e && s <= end;
		return overlap;
	}


	public boolean intersect(SAMRecord sr) {

		int s = sr.getAlignmentStart();
		int e = sr.getAlignmentEnd();
		if(e==0||e==s){
			e = sr.getAlignmentStart() + sr.getReadLength();
		}
		boolean overlap = (start <= e && s <= end);
		return overlap;

	}
	
	public boolean equals(Object obj){
		
		if(!(obj instanceof CapInterval)){
			return false;
		}
		CapInterval ci0 = (CapInterval)obj;
		return ci0.chr.equals(chr) 
			&& (ci0.start == start)
			 && (ci0.end == end);
	}
	
	public int hashCode(){
		return this.infoStr().hashCode();
	}

	public String infoStr() {
		return chr+":"+start+"-"+end+" ("+(length)+")";
	}

	CNVInfo cnvinfo;
	public void setCNVInfo(CNVInfo _cnvinfo) {
		cnvinfo = _cnvinfo;
	}
	public CNVInfo getCNVInfo(){
		return cnvinfo;
	}

	public void setDenioseValue(double denoise) {
		cnvinfo.setDenioseValue(denoise);
	}

	public void setCN(double copynumber) {
		cnvinfo.setCN(copynumber);		
	}

	public double getOriginalValue() {
		return cnvinfo.getOriginalTnratio();
	}
	
	public double getValue() {
		return cnvinfo.getTnratio();
	}

	public boolean include(int pos) {
		boolean overlap = (start <= pos && pos <= end);
		return overlap;
	}

	public boolean include(int pos, int mergin) {
		boolean overlap = (start-mergin <= pos && pos <= end+mergin);
		return overlap;
	}

	
	public double getDenioseValue() {
		return cnvinfo.getDenoise();
	}

	public double getCN() {
		return cnvinfo.getCopynumber();
	}

	double hmmvalue=0;
	double varidateval = 0;
	public void setHMMValue(double d) {
		 hmmvalue=d;
		 varidateval=d;
	}

	public double getHMMValue() {
		return  hmmvalue;
	}

	public float getDuality() {
		return duality;
	}

	public void setGCAdjustedTNratio(double adjustedY) {
		if(cnvinfo!=null){
			cnvinfo.tnratio = adjustedY;
		}
	}

	public double getVaridateVal() {
		
		return varidateval;
	}

	public void setVaridateval(double varidateval) {
		this.varidateval = varidateval;
	}



	
	

}