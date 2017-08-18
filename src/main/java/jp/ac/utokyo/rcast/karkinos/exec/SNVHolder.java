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

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPBean;
import jp.ac.utokyo.rcast.karkinos.annotation.stats.Fisher;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.graph.output.FormatHelper;

public class SNVHolder implements java.io.Serializable{
	

	String chr;
	int pos;
	CapInterval ci;
	DbSNPBean dbSNPbean;
	FilterResult filterResult;	
	
	double pvalFisher = 0;
	boolean fisherTestSignif = false;		
	
	public double getPvalFisher() {
		return pvalFisher;
	}

	public void setPvalFisher(double pvalFisher) {
		this.pvalFisher = pvalFisher;
	}

	public boolean isFisherTestSignif() {
		return fisherTestSignif;
	}

	public void setFisherTestSignif(boolean fisherTestSignif) {
		this.fisherTestSignif = fisherTestSignif;
	}

	public FilterResult getFilterResult() {
		return filterResult;
	}

	public void setFilterResult(FilterResult filterResult) {
		this.filterResult = filterResult;
		filterResult.setDbSNPbean(dbSNPbean);
	}

	public String getInfoStr() {
		
		StringBuffer sb = new StringBuffer();
		sb.append(chr+"\t");
		sb.append(pos+"\t");
		String dbSNP = "";
		if(dbSNPbean!=null){
			dbSNP = dbSNPbean.getInfo();
		}
		sb.append(dbSNP+"\t");
		sb.append(ci.getCN()+"\t");
		sb.append(normal.getInfoStr()+"\t");
		sb.append(tumor.getInfoStr()+"\t");
		
		
		return sb.toString();
	}
	
	public DbSNPBean getDbSNPbean() {
		return dbSNPbean;
	}
	public void setDbSNPbean(DbSNPBean dbSNPbean) {
		this.dbSNPbean = dbSNPbean;
	}
	public CapInterval getCi() {
		if(ci==null){
			return new CapInterval("chr0", 0, 0, false);
		}
		return ci;
	}
	public void setCi(CapInterval ci) {
		this.ci = ci;
	}
	public int getPos() {
		return pos;
	}
	public void setPos(int pos) {
		this.pos = pos;
	}
	public String getChr() {
		return chr;
	}
	int flg;
	PileUPResult normal;
	PileUPResult tumor;
	
	public int getFlg() {
		return flg;
	}
	public void setFlg(int flg) {
		this.flg = flg;
	}
	public PileUPResult getNormal() {
		return normal;
	}
	public void setNormal(PileUPResult normal) {
		this.normal = normal;
	}
	public PileUPResult getTumor() {
		return tumor;
	}
	public void setTumor(PileUPResult tumor) {
		this.tumor = tumor;
	}
	public void setChr(String _chr) {
		chr = _chr;		
	}

	boolean hetroSNP = false;
	public void serHetroSNP(boolean _hetroSNP) {
		hetroSNP = _hetroSNP;	
	}

	public boolean isHetroSNP() {
		return hetroSNP;
	}

	
	

	
}
