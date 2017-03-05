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

import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class CNVUtils {

	public static void reflectToSNV(DataSet dataset,
			List<CopyNumberInterval> lohs, List<CopyNumberInterval> gains) {

		List<CopyNumberInterval> hdAmp = new ArrayList<CopyNumberInterval>();
		int unit = 5000000;

		for (CopyNumberInterval cni : dataset.getCniVaridateList()) {

			if (cni.getCopynumber() == 2)
				continue;
			//
			// int length = Math.abs(cni.getEnd()-cni.getStart());
			// boolean amp = cni.getCopynumber()>2;
			// boolean narrowamp = amp&&(length<unit);
			// if(cni.isAllelic()&&cni.getNoSNP()>100&&(length>unit)&&!narrowamp){
			//
			//
			// al.add(cni);
			// }
			if (cni.getCopynumber() == 0 ) {
				hdAmp.add(cni);
			}
			if(focalamp(cni)){
				hdAmp.add(cni);
			}
		}

		List<List<WaveletIF>> cap = dataset.getCapInterval();
		for (List<WaveletIF> list : cap) {

			for (WaveletIF wi : list) {
				CapInterval ci = (CapInterval) wi;
				CopyNumberInterval ic1 = intercect(lohs, ci, false);
				// set UPD
				if (ic1 != null) {

					CopyNumberInterval icgain = intercect(gains, ci, true);
					if (icgain != null) {

						if (ci.getCnvtotal() == 2.0) {
							float sumallele = ci.getAafreq() + ci.getBafreq();
							ci.setBafreq(0);
							ci.setAafreq(sumallele);
						}
					}
				}

				CopyNumberInterval ic2 = intercect(hdAmp, ci, true);
				if (ic2 != null) {
					ci.setVaridateval(ic2.getCopynumber() / 2);
				}

			}
		}

	}

	private static boolean focalamp(CopyNumberInterval cni) {
		
		boolean b1 = cni.getCopynumber() >=6;
		boolean b2 = Math.abs(cni.getEnd()-cni.getStart()) < 1000000;
		return b1 && b2;
		
	}

	private static CopyNumberInterval intercect(List<CopyNumberInterval> lohs,
			CapInterval ci, boolean gain) {

		CopyNumberInterval ret = null;
		for (CopyNumberInterval cil : lohs) {

			if (gain) {
				if (cil.getCopynumber() < 1)
					continue;
			} else {
				if (cil.getCopynumber() > 1)
					continue;
			}
			if (!cil.getChr().equals(ci.getChr())) {
				continue;
			}
			if (cil.getStart() <= ci.getEnd() && ci.getStart() <= cil.getEnd()) {
				return cil;
			}

		}
		return ret;

	}

}
