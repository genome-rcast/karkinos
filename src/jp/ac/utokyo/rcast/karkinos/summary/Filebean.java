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

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Filebean {

	File tdata = null;
	File csvdata = null;
	File pdfdata = null;
	
	public File getPdfdata() {
		return pdfdata;
	}
	public void setPdfdata(File pdfdata) {
		this.pdfdata = pdfdata;
	}
	String id = "";
	
	
	public List<TCBean> tcbeanL = new ArrayList<TCBean>();
	int tcflg = 0;
	public int getTcflg() {
		return tcflg;
	}
	public void setTcflg(int tcflg) {
		this.tcflg = tcflg;
	}
	public void addTC(String name, float tr, float sd, int nohsnp, float correl) {
		//
		TCBean bean = new TCBean();
		bean.setName(name);
		bean.setTumorratio(tr);
		bean.setNosnp(nohsnp);
		bean.setSd(sd);
		bean.setCorrel(correl);
		
		tcbeanL.add(bean);
	}
	
	public TCBean getTCBean() {
		if(tcbeanL.size()==0){
			return null;
		}
		if(tcflg>0&&(tcbeanL.size()>=tcflg)){
			return tcbeanL.get(tcflg-1);
		}else{
			return tcbeanL.get(0);
		}
	}
	
	public List<CNAInterval> cnalist;
	public Map<String, String[]> datamap;
	
	
	
	public File getTdata() {
		return tdata;
	}
	public void setTdata(File tdata) {
		this.tdata = tdata;
	}
	public File getCsvdata() {
		return csvdata;
	}
	public void setCsvdata(File csvdata) {
		this.csvdata = csvdata;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	float tr = 0f;
	public void settr(float tr) {
			this.tr = tr;	
	}
	public float getTr() {
		return tr;
	}
	
}
