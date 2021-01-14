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

import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class GeneEachCNV {
	public static void calcCNVForEachgene(DataSet dataset, GeneExons ge) {
		List<List<WaveletIF>> cap = dataset.getCapInterval();

		long normaltotal=0;
		long tumortotal=0;

		for (List<WaveletIF> list : cap) {
			for (WaveletIF wi : list) {
				CapInterval civ = (CapInterval) wi;
				String chr = civ.getChr();
				int start = civ.getStart();
				int end = civ.getEnd();

				String gen = ge.getGeneId(chr, (int)start);

				if(gen==null){
					gen = ge.getGeneId(chr, (int)end);
				}

				long normalcnt = civ.getCNVInfo().getNormalcnt();
				long tumorcnt = civ.getCNVInfo().getTumorcnt();

				if(gen!=null){
				 normaltotal = normaltotal+normalcnt;
				 tumortotal = tumortotal+tumorcnt;
				}

				if(gen!=null){
					DataHolder dh = ge.getCounterForGene().get(gen);

					if(dh==null){
						dh = new DataHolder();
						ge.getCounterForGene().put(gen, dh);
					}
					if(normalcnt > 50){
							dh.put((int)normalcnt,(int)tumorcnt,civ );
					}
				}
			}
			ge.setNormaltotal(normaltotal);
			ge.setTumortotal(tumortotal);
		}
	}
}
