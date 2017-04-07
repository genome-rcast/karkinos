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

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Float;
import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;

public class CorrelVaridate {

	final static int FINALV =9;
	public static void recheck(DataSet dataset) {
		
		List<CopyNumberInterval> cniListb4 = dataset.getCniVaridateList();
		ArrayList<CopyNumberInterval> cniListb4c = new ArrayList<CopyNumberInterval>();
		cniListb4c.addAll(cniListb4);
		
		List<CopyNumberInterval> cniList = dataset
		.getCopyNumberIntervalList(FINALV);
		
		//HDcheck
		CopyNumberInterval cnib4 =null;
		for(CopyNumberInterval cni: cniList){
			
			//
			if(!cni.isHdelation()){
				if(Math.abs(cni.getEnd()-cni.getStart())<1000000){
					if(cni.getCopynumber()>6){
						
						if(cnib4!=null){
							if((cnib4.getCopynumber()-cni.getCopynumber())>1){
								cni.setHdelation(true);
							}
						}
						
					}
				}
			}
			cnib4 = cni;
			
		}
		
		dataset.setVaridatedCNlist(cniList);
		
		for(CopyNumberInterval cni:cniListb4c){
			if(cni.isHdelation()){
				if(!cniList.contains(cni)){
					cniList.add(cni);
				}
			}			
		}
		
	}
	
	//
	public static void varidate(DataSet dataset, AllelicCNV alCNV) {

		List<CopyNumberInterval> cniList = dataset
				.getCopyNumberIntervalList(DataSet.MODE_HMM);
		//to integrate allelic and nonallelic CNA call
		//
		List<CopyNumberInterval> allelicCNAs = alCNV.getSortedList();
		
		
		for (CopyNumberInterval cni : cniList) {

			List<Point2D.Float> list = new ArrayList<Point2D.Float>();
			List<SNVHolder> holderL = new ArrayList<SNVHolder>();
			
			if(include(allelicCNAs ,cni)){
				cni.setSupportbyAllelic(true);
			}
			
			//
			int noSNP = 0;
			for (SNVHolder holder : dataset.getSnvlist()) {

				if (!include(holder, cni))
					continue;
				
				holderL.add(holder);
				if (holder.isHetroSNP()) {

					float nr = holder.getNormal().getRatio();
					float tr = holder.getTumor().getRatio();
					list.add(new Point2D.Float(nr, tr));
					noSNP++;					
				}

			}
			cni.setNoSNP(noSNP);
			float snpcorrel = (float)getPearsonCorrelation(list);
//			System.out.println("start");
//			for(Point2D d:list){
//				System.out.println(d.getX()+"\t"+d.getY());
//			}
			//System.out.println("correl="+snpcorrel);
			
			cni.setSnpclrrel(snpcorrel);
			//correlation corpse because of CNA
			boolean varidated = snpcorrel < 0.7;
			if(noSNP<50){
				varidated = true;// cannot varidate
			}
			cni.setVaridated(varidated);

			//
			for(SNVHolder holder :holderL){
				//
				if(!varidated){
					double hmmv = holder.getCi().getHMMValue();
					if(1<hmmv && hmmv<=1.5){
						holder.getCi().setVaridateval(1);
					}
					if(1>hmmv && hmmv>=0.5){
						holder.getCi().setVaridateval(1);
					}
				}
			}
			
		}

		//check posiible amplification/delation
				
		
		
		//
//		for(CopyNumberInterval cni:lowList){
//			//
//			if(!cni.isRecurrent()){
//				cniList.add(cni);
//			}
//		}
		
		cniList.addAll(allelicCNAs);
		// set to dlist
		dataset.setVaridatedCNlist(cniList);

	}



	private static boolean include(List<CopyNumberInterval> allelicList,
			CopyNumberInterval cni) {
		for(CopyNumberInterval allelic:allelicList){
			
			if(include(cni,allelic)){
				return true;
			}
		}		
		return false;
	}



	private static boolean include(CopyNumberInterval cni, CopyNumberInterval allelic) {
		if(!cni.getChr().equals(allelic.getChr())){
			return false;
		}
		
		if((cni.getStart()<allelic.getEnd())&&(allelic.getStart()<cni.getEnd())){
			
			if(allelic.getCopynumber()<2){
				if(cni.getCopynumber()<2){
					 cni.setRecurrent(true);
					 return true;
				}
			}
			
			if(allelic.getCopynumber()>2){
				if(cni.getCopynumber()>2){
					 cni.setRecurrent(true);
					 return true;
				}
			}
			
		}		
		
		
		return false;
	}



	public static double getPearsonCorrelation(List<Float> list) {
		
		if(list.size()==0)return 0;
		double result = 0;
		double sum_sq_x = 0;
		double sum_sq_y = 0;
		double sum_coproduct = 0;
		double mean_x = list.get(0).getX();
		double mean_y = list.get(0).getY();
		for (int i = 2; i < list.size() + 1; i += 1) {
			double sweep = (double) (i - 1) / i;
			double delta_x = list.get(i-1).getX() - mean_x;
			double delta_y = list.get(i-1).getY() - mean_y;
			sum_sq_x += delta_x * delta_x * sweep;
			sum_sq_y += delta_y * delta_y * sweep;
			sum_coproduct += delta_x * delta_y * sweep;
			mean_x += delta_x / i;
			mean_y += delta_y / i;
		}
		double pop_sd_x = Math.sqrt(sum_sq_x / list.size());
		double pop_sd_y = Math.sqrt(sum_sq_y / list.size());
		double cov_x_y = sum_coproduct /list.size();
		result = cov_x_y / (pop_sd_x * pop_sd_y);
		return result;
	}

	private static boolean include(SNVHolder holder, CopyNumberInterval cni) {

		if (holder.getChr().equals(cni.getChr())) {
			if (cni.getStart() <= holder.getPos()
					&& cni.getEnd() >= holder.getPos()) {
				return true;
			}
		}
		return false;
	}



	

}
