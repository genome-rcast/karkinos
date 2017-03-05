package jp.ac.utokyo.rcast.karkinos.utils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class GeneEachCNV {

	public static void calcCNVForEachgene(DataSet dataset, GeneExons ge) {
		
		
		List<List<WaveletIF>> cap = dataset.getCapInterval();
		//
		
		
		long normaltotal=0;
		long tumortotal=0;
		
		for (List<WaveletIF> list : cap) {

			for (WaveletIF wi : list) {

				CapInterval civ = (CapInterval) wi;
				String chr = civ.getChr();
				int start = civ.getStart();
				int end = civ.getEnd();
				
				//
				String gen = ge.getGeneId(chr, (int)start);
				//
				if(gen==null){
					gen = ge.getGeneId(chr, (int)end);
				}
				//				
				//
				int length = civ.getLength();
				long normalcnt = civ.getCNVInfo().getNormalcnt();
				long tumorcnt = civ.getCNVInfo().getTumorcnt();
				
				if(gen!=null){
				 normaltotal = normaltotal+normalcnt;
				 tumortotal = tumortotal+tumorcnt;
				}
				
				if(gen!=null){
					
					//
					DataHolder dh = ge.getCounterForGene().get(gen);
					//
					if(dh==null){
						//
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
