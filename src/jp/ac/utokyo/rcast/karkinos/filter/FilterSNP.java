package jp.ac.utokyo.rcast.karkinos.filter;

import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;

public class FilterSNP {

	public static String filter(PileUPResult pir) {
		
		//
		//qual by depth
		double ph = (pir.getPhred()/10);
		int altcnt = pir.getAltCnt();
		
		double QD = (double)ph/altcnt;
		if(QD<2){
			return "LowQD";
		}
		
		if(pir.getRatio()<0.2){
			return "LowAF";
		}
		
		if(pir.getMQrms()<40){
			return "LowMQ";
		}		
		//FS
		
			
		//MQ
		
		//
		
		
		return "PASS";
	}

	public static String filterIndel(PileUPResult pir) {
		
		if(pir.getRatio()<0.2){
			return "LowAF";
		}
		
		if(pir.getMQrms()<40){
			return "LowMQ";
		}
		
		return "PASS";
	}

	

}
