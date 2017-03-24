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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import jp.ac.utokyo.rcast.karkinos.filter.BaseandQData;

public class OddRatioCalculator {

	public static float getAdjuatedLogt(char genomeR, char alt,
			List<BaseandQData> pileuplist, double tumorratio, double copynumber,double normalAF) {

		//int copynumbern = getRound(copynumber);
		
		List<BaseandQData> pileuplistToCheck = getTocheckList(genomeR, alt,
				pileuplist, tumorratio, copynumber,normalAF);
		// /
		double reflikehood = 0;
		double altlikehood = 0;
		for (BaseandQData data : pileuplistToCheck) {

			double qual0 = (int) data.getQual() & 0xFF;
			qual0 = qual0 * 0.1;
			double pNomatch = (1 / Math.pow(10, qual0));
			double pmathch = 1 - pNomatch;

			if (equalChar(data.getBase(), genomeR)) {

				reflikehood = reflikehood + pmathch;
				altlikehood = altlikehood + (pNomatch / (double) 3);

			} else {

				reflikehood = reflikehood + (pNomatch / (double) 3);
				altlikehood = altlikehood + pmathch;

			}

		}
		double d = altlikehood / reflikehood;
		return (float) (Math.log10(d));

	}


	private static List<BaseandQData> getTocheckList(char genomeR, char alt,
			List<BaseandQData> pileuplist, double tumorratio, double copynumber, double normalAF) {

		//

		List<BaseandQData> reflist = new ArrayList<BaseandQData>();
		List<BaseandQData> altlist = new ArrayList<BaseandQData>();
		//
		for (BaseandQData data : pileuplist) {
			//
			if (equalChar(data.getBase(), genomeR)) {
				reflist.add(data);
			} else if (equalChar(data.getBase(), alt)) {
				altlist.add(data);
			}

		}

//		int addnum = getAddSampleNumber(pileuplist.size(), altlist.size(),
//				tumorratio, copynumber);
//		int cnt = 0;
//		for (BaseandQData data : refset) {
//			if (cnt > addnum - 1) {
//				break;
//			}
//			//
//			altlist.add(data);
//			cnt++;
//		}
		List<BaseandQData> retlist = new ArrayList<BaseandQData>();
		int cntref = (int)(reflist.size()*tumorratio);
		int cntalt = altlist.size() - (int)Math.floor((pileuplist.size()*normalAF));
		for(int n=0;n<cntref;n++){
			retlist.add(reflist.get(n));
		}
		for(int n=0;n<cntalt;n++){
			retlist.add(altlist.get(n));
		}
		return retlist;
	}

	static final float noize_factor = 0.8f;
	private static int getAddSampleNumber(int totaldepth, int altsize,
			double tumorratio, double copynumber) {
		
		double adjusttotal = (totaldepth * tumorratio * noize_factor);
		int sampleExpectNum = (int) (adjusttotal / (double) copynumber);
		if (altsize < sampleExpectNum) {
			return sampleExpectNum - altsize;
		}
		return 0;
		
	}

	private static boolean equalChar(char c1, char c2) {

		c1 = Character.toUpperCase(c1);
		c2 = Character.toUpperCase(c2);
		return c1 == c2;
	}

}
