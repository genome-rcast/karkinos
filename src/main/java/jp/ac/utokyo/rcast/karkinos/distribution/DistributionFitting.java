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
package jp.ac.utokyo.rcast.karkinos.distribution;

import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class DistributionFitting {

	public static float[] getObserveRatio(DataHolderByCN dataHolderByCN,
			float degree, double tumorvariance) {
		if (dataHolderByCN == null) {
			return new float[]{0,0};
		}
		if(tumorvariance>5){
			tumorvariance=5;
		}
		//
		int maxratio = 0;
		double maxcp = 0;
		for (int n = 50; n <= 100; n++) {
			//
			float[] theoreticalDist = getThereticalDist(n, degree,
					tumorvariance);
				
			double cp = crossProduct(theoreticalDist, dataHolderByCN.snpDistT);
			System.out.println("n="+n+"cp="+cp);
			if (cp > maxcp) {
				maxratio = n;
				maxcp = cp;
			}

		}
		return new float[]{maxratio,(float) maxcp};
	}
	
	public static float[] getObserveRatio(DataHolderByCN dataHolderByCN,
			float degree, double tumorvariance,int start,int end) {
		if (dataHolderByCN == null) {
			return new float[]{0,0};
		}

		//
		if(tumorvariance>5){
			tumorvariance=5;
		}
		
		int maxratio = 0;
		double maxcp = 0;
		for (int n = start; n <= end; n++) {
			//
			float[] theoreticalDist = getThereticalDist(n, degree,
					tumorvariance);
				
			double cp = crossProduct(theoreticalDist, dataHolderByCN.snpDistT);
			System.out.println("n="+n+"cp="+cp);
			if (cp > maxcp) {
				maxratio = n;
				maxcp = cp;
			}

		}
		return new float[]{maxratio,(float) maxcp};
	}
	

	private static double crossProduct(float[] theoreticalDist, float[] snpDistT) {

		// exclude < 2 and > 98 : both edge
		// exclude 47 to 53 : 6% center 
		float[] normalizesnpDistT = norm(snpDistT);
		double d = 0;
		for (int n = 2; n < 98; n++) {
			//
			if(n>=47 && n<=53){
				continue;
			}
			d = d + (theoreticalDist[n] * normalizesnpDistT[n]);
		}
		return d;
	}

	private static float[] norm(float[] dist) {
		float sum = 0;
		int idx = 0;
		for(float f:dist){
			if(idx==0){
				idx++;
				continue;
			}
			if(idx==99){
				idx++;
				continue;
			}
			if(idx>=45&&idx<=55){
				idx++;
				continue;
			}
			sum = sum+f;
			idx++;
		}
		float[] ndist = new float[dist.length];
		
		//
		int n= 0;
		for(float f:dist){
			if(n==0){
				n++;
				continue;
			}
			if(n==99){
				n++;
				continue;
			}
			if(n>=45&&n<=55){
				n++;
				continue;
			}
			ndist[n] = (float)((double)f/(double)sum);
			n++;
		}				
		return ndist;
	}

	private static float[] getThereticalDist(int ratio, float degree,
			double variance) {

		float x = 0;
		double sum = 0;
		float[] dist = new float[100];
		for (int n = 0; n < 100; n++) {

			x = n + 0.5f;
			float y = dist(x, variance, ratio);
			sum = sum+y;
			dist[n]=y;
		}
		int idx =0;
		for(double d:dist){
			
			dist[idx] = (float) (dist[idx]/sum);
			idx++;
		}
		return dist;
	}

	private static float dist(float x, double sd, int ratio) {
		
		if(sd==0){
			sd = 10;// to avoid error
		}

		if (ratio == 50) {
			float mean = 50;
			NormalDistribution normald = new NormalDistributionImpl(mean, sd);
			return (float) normald.density((double) x);
		} else {

			NormalDistribution normald1 = new NormalDistributionImpl(ratio, sd);
			NormalDistribution normald2 = new NormalDistributionImpl(
					50 - (ratio - 50), sd);
			return (float) (normald1.density((double) x) + normald2
					.density((double) x));

		}

	}

}
