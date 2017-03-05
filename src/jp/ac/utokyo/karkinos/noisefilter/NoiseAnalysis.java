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
package jp.ac.utokyo.karkinos.noisefilter;

import static jp.ac.utokyo.rcast.karkinos.filter.FilterResult.INFO_SUPPORTED_BY_ONEDirection;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.FilterResult;
import jp.ac.utokyo.rcast.karkinos.filter.SupportReadsCheck;
import jp.ac.utokyo.rcast.karkinos.readssummary.SNPDepthCounter;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;
import jp.ac.utokyo.rcast.karkinos.utils.GenotypeKeyUtils;

public class NoiseAnalysis {
	//

	AFDepthMatrix afm = null;
	List<BinData> binlist = null;
	
	int snvcount=0;
	int atocCount = 0;
	int lowcount =0;
	boolean highErrorSample = false;
	public void analysisNoiseRegion(DataSet dataset, int ploidy) {

		int binsize = 20;
		
		snvcount=0;
		atocCount=0;
		lowcount=0;

		float tumorContentsRatio = dataset.getAnalyseDist().getTumorratio();
		if(tumorContentsRatio == 0){
			tumorContentsRatio = dataset.getTumorRatio();
		}

		binlist = getBorderDepth(dataset.getSnvlist(), tumorContentsRatio,
				binsize);
 
		// from hetroSNP
		afm = new AFDepthMatrix(binlist, binsize);
		for (SNVHolder snv : dataset.getSnvlist()) {

			float f = (float) snv.getCi().getVaridateVal();
			//
			// if (ploidy == 4) {
			// if (!(f == 2.0)) {
			// continue;
			// }
			// }else{
			// if (!(f == 1.0)) {
			// continue;
			// } 
			// }

			if (!((f == 1.0)||(f == 2.0))) {
				continue;
			}
			if (snv.isHetroSNP()) {

				afm.regHetroSNP(snv, tumorContentsRatio);

			}

		}

		for (SNVHolder snv : dataset.getSnvlist()) {

			float f = (float) snv.getCi().getVaridateVal();
			if (!((f == 1.0)||(f == 2.0))) {
				continue;
			}
			int flg = snv.getFlg();
			boolean check = (flg == PileUP.SomaticMutation || snv.getFlg() == PileUP.TumorINDEL);
			if (check && snv.getFilterResult() != null
					&& (snv.getFilterResult().isPassFilter() 
							||(snv.getFilterResult().getPassFilterFlg().contains(FilterResult.illuminaSpecific))) ) {

				afm.regSNV(snv, tumorContentsRatio);
				snvcount++;
				char genomeR = snv.getTumor().getGenomeR();
				char altTumor = snv.getTumor().getALT();
				boolean lowsuport = snv.getTumor().getAltCnt()<=4;
				if(lowsuport){
					lowcount++;
				}
				
				if (genomeR == 'C' && altTumor == 'A') {
					atocCount++;
				}

				if (genomeR == 'G' && altTumor == 'T') {
					atocCount++;					
				}
				
			}

		}
		// sort snv list for SNV cal
		afm.sortList();
		// excute em
		afm.EMmethod(tumorContentsRatio);
		// find regression
		afm.calcregresion();
		double atoCratio = ((double)atocCount)/((double)snvcount);
		double lowratio = ((double)lowcount)/((double)snvcount);
		
		highErrorSample = (snvcount>=500)&& ((atoCratio>0.85)||(lowratio>0.8));
		
	}

	public List<BinData> getBinlist() {
		return binlist;
	}

	public AFDepthMatrix getAfm() {
		return afm;
	}

	private List<BinData> getBorderDepth(List<SNVHolder> snvlist,
			float tumorContentsRatio, int binsize) {

		List<BinData> binlist = new ArrayList<BinData>();
		Map<Integer, NCounter> counter = new HashMap<Integer, NCounter>();
		int total = 0;
		int depthmax = 0;
		for (SNVHolder snv : snvlist) {
			float f = (float) snv.getCi().getVaridateVal();
			if (!((f == 1.0)||(f == 2.0))) {
				continue;
			}
			if (snv.isHetroSNP()) {

				total++;
				int depth = snv.getTumor().getTotalcnt();
				if (depth > depthmax) {
					depthmax = depth;
				}
				if (counter.containsKey(depth)) {
					counter.get(depth).inc();
				} else {
					NCounter nc = new NCounter();
					counter.put(depth, nc);
				}

			}
		}
		int unit = total / binsize;
		int cnt = 0;
		int binstart = 0;
		for (int n = 1; n < depthmax; n++) {

			if (counter.containsKey(n)) {

				if (binstart == 0) {
					binstart = n;
				}
				//
				cnt = cnt + counter.get(n).getN();
				if (cnt > unit) {
					binlist.add(new BinData(binstart, n));
					cnt = 0;
					binstart = 0;
				}

			}

		}
		return binlist;

	}

	public boolean reject(Point2D p2d) {

		return afm.reject((int) p2d.getY(), (float) p2d.getX(),highErrorSample);
	}

	public float getPval(float adjustedratio, float readdepth) {

		return afm.getPval(adjustedratio, readdepth);

	}

	public boolean reject(SNVHolder snv, float tumorContentsRatio, boolean indel) {

		float adjusttedAF = CalcUtils.getTumorrateAdjustedRatio(snv,
				tumorContentsRatio);
		int depth = snv.getTumor().getTotalcnt();
		int depth_c = (int) (depth * tumorContentsRatio);
		
		boolean oxoGCand = SupportReadsCheck.oxoG(snv.getTumor().getGenomeR(), snv.getTumor().getALT());
		boolean ffpeCand = SupportReadsCheck.ffpe(snv.getTumor().getGenomeR(), snv.getTumor().getALT());
		boolean onedirec = snv.getFilterResult().getInfoFlg().contains(INFO_SUPPORTED_BY_ONEDirection);
		
		
		if((oxoGCand||ffpeCand ) && onedirec){
			return afm._reject(depth_c, adjusttedAF);
		}
		
		return afm.reject(depth_c, adjusttedAF,highErrorSample);
		
		

		// if(indel){
		// //since indel is hard to detect
		// depth_c = depth_c*2;
		// }
		// char ref = snv.getNormal().getGenomeR();
		// char alt = snv.getNormal().getALT();
		// String key = ref + "to" + alt;
		// key = GenotypeKeyUtils.aggrigateKeys(key);
		// if(key.equals("")){
		// key = "Indel";
		// }
		// AFDepthMatrix afm = null;
		// if(matrix.containsKey(key)){
		// afm = matrix.get(key);
		// return afm.reject(depth_c,adjusttedAF);
		// }
		// return AFDepthMatrix.defultreject(depth_c, adjusttedAF);
	}

}
