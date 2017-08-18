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
package jp.ac.utokyo.karkinos.ploidy;

public class VirtualPeakCalculator {

	public static final int diploid = 2;
	public static final int tetraploid = 4;

	public static void calc(int[] ploidy, int baseploidy, float[] imbalance,
			float[] distToBase,double peakpos) {

		for (int n = 0; n < 100; n++) {

			int tumorpurity = n;
			float disttobase = getDistTobase(ploidy, baseploidy, tumorpurity, peakpos);
			float imblance = getImbalance(ploidy, baseploidy, tumorpurity);
			imbalance[n] = imblance;
			distToBase[n] = disttobase;

		}

	}

	private static float getImbalance(int[] ploidy, int baseploidy,
			int tumorpurity) {

		double pt = 0.01*tumorpurity;
		//even peak
		if (Math.abs(ploidy[0]-ploidy[1]) == 0) {
			
				return 0f;
			
		} else if (ploidy[1] >  ploidy[0]) {

				//shoul not happen
			return 0f;
				
		} else  {

			//r = (m-1)x+1 / ((m+l)-2)x +2
			int c1 =  ploidy[0]-1;
			int c2 = ploidy[0]+ploidy[1]-2;
			double part1 = (c1*pt)+1;
			double part2 = (c2*pt)+2;
			double r = part1/part2;
			float im = (float)(r-0.5);			
			if(im<0){
				im=0f;
			}
			return im;
		}

	}

	private static float getDistTobase(int[] ploidy, int baseploidy,
			int tumorpurity,double peakpos) {

		int totalpolidy = ploidy[0]+ploidy[1];
		double pt = 0.01*tumorpurity;
		float unit = getUnitDist(pt,baseploidy,peakpos);
		if (baseploidy == diploid) {

			if(totalpolidy==diploid){
				return 0;
			}else{
				
				int diff = totalpolidy-diploid;
				return unit*diff;
				
			}
			
			
		} else {
			//tetraploid
			if(totalpolidy==tetraploid){
				return 0;
			}else{
				
				int diff = totalpolidy-tetraploid;
				return (unit*diff);
				
			}
			
		}

	}

	private static float getUnitDist(double pt, int baseploidy, double peakpos) {
		
		if (baseploidy == diploid) {
			return (float)((pt/2)*peakpos);
		}else{
			//y = -0.15x2+0.4x
			return (float) ((-0.15*Math.pow(pt, 2) + (0.4*pt)) * peakpos) ;
		}
		
	}

}
