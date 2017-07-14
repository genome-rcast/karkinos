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
package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.util.ArrayList;
import java.util.List;

public class PeakDistResolve {
	

	public static List<Double> resolvePeakDist(PeaksInfo pi) {
		
		//find peak dist
		List<Double> peakdist = new ArrayList<Double>();
		List<Peak> peaklist = pi.getCopy();
		
		
		
		
		for (int m= 0;m<peaklist.size()-1;m++) {
			
			Peak p1 = peaklist.get(m);
			Peak p2 = peaklist.get(m+1);
			double dist = Math.abs(p1.getU()-p2.getU());
			peakdist.add(dist);
						
		}		
		
		List<Double> peakdistcand = new ArrayList<Double>();
		//
		
		
		//
		return peakdistcand;
		
		
	}

}
