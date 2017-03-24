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
package jp.ac.utokyo.rcast.karkinos.utils;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;

public class CalcUtils {
	
	public static float getTumorrateAdjustedRatio(SNVHolder snv,float tratio) {

		
		
		double nontumorr = 1 - tratio;
		double nAF = snv.getNormal().getRatio();
		double tAF = snv.getTumor().getRatio();
		
		
		
		
//		if (nontumorr > 1 || nontumorr < 0) {
//			// something wrong
//			return (float) tAF;
//		} else {
//
//			float f = (float) (tAF - (1-tratio)*nAF)/tratio;
//			if(f>1){
//				f=1;
//			}
//			return f;
//			
//		}
		//float aaf = (float) (snv.getCi().getAafreq());
		//float baf = (float) (snv.getCi().getBafreq());
		float taf = (float) (snv.getCi().getCnvtotal());
		//
		if (nontumorr > 1 || nontumorr < 0) {
			// something wrong
			return (float) tAF;
		} else {

			float f = (float) ((tAF - nAF)*((tratio*taf)+(2*nontumorr)))/(tratio*taf);
			if (f > 1) {
				f = 1;
			}
			return f;

		}
		
		
	}
	
	public static String revcon(String read) {

		StringBuffer sb = new StringBuffer();
		for (char c : read.toCharArray()) {
			sb.append(comp(c));
		}
		sb.reverse();
		return sb.toString();
	}

	private static char comp(char c) {

		if (c == 'A') {
			return 'T';
		} else if (c == 'T') {
			return 'A';
		} else if (c == 'C') {
			return 'G';
		} else if (c == 'G') {
			return 'C';
		}
		return 'N';
	}

	
	public static boolean pass2(FilterResult filterResult) {
		
		if(!filterResult.isPassFilter()){
			return false;
		}else{
			//cosmic
			for(int flg:filterResult.getInfoFlg()){
			  if(flg==FilterResult.INFO_COSMIC_Validate){
				return true;
			  }	
			}			
			// if Indel
			
		}
		for(int flg:filterResult.getInfoFlg()){
			if(flg>200){
				return false;
			}
		}		
		return true;
	
	}
}
