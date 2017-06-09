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
package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.wavelet.FunctionRegression;

public class SNVHolderPlusACnv {
	
	
	public SNVHolderPlusACnv(SNVHolder snv, double ntratio){
		this.ntratio = ntratio;
		this.snv = snv;
		setVal();
	}
	
	public SNVHolder getSnv() {
		return snv;
	}

	

	boolean valid = false;
	double ntratio = 0;
	
	private void setVal() {
		
		//
		int refidx = PileUPResult.seqALL.indexOf(snv.getNormal().getGenomeR());
		if(refidx<0){
			return;
		}
		valid = true;
		int normalref = snv.getNormal().getSeqCounter()[refidx]+1;
		int normalalt = snv.getNormal().getAltCnt()+1;
		
		char alt = snv.getNormal().getALT();
		int altidx = PileUPResult.seqALL.indexOf(alt);
		int tumorref = snv.getTumor().getSeqCounter()[refidx]+1;
		int tumoralt = snv.getTumor().getSeqCounter()[altidx]+1;
		
		int normaltotal = snv.getNormal().getTotalcnt();
		int tumortotal = snv.getTumor().getTotalcnt();
		
//		float tr = ratio(tumortotal,normaltotal);		
//		if(tr>1.25||tr>0.75){
//			tr = 1;
//		}
		float r1 = ratio(tumorref,normalref,ntratio);
		float r2 = ratio(tumoralt,normalalt,ntratio);
		
		lowera = new ACNVInfoBean();
		highera = new ACNVInfoBean();
		
		if(r1<r2){
			lowera.row = r1;
			highera.row = r2;
		}else{
			lowera.row = r2;
			highera.row = r1;
		}
		
		
	}

	private float ratio(double a, double b, double ntratio2) {
		
		float r = (float)((a/b)/ntratio2);
		if(r>6){
			r=6;
		}
		return r;
	}	

	public ACNVInfoBean getLowera() {
		return lowera;
	}

	public ACNVInfoBean getHighera() {
		return highera;
	}



	FunctionRegression fr;
	SNVHolder snv;
	ACNVInfoBean lowera;
	ACNVInfoBean highera;
	

}
