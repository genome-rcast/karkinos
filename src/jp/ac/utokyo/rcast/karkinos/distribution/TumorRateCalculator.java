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

import java.util.Map;

public class TumorRateCalculator {

	private static final float LOH = 0.5f;
	private static final float GAIN = 1.5f;
	private static final float MINForGAIN = 0.33f;
	private static final float MAXForGAIN = 0.66f;


	public static float getTumorRatio(float degree, float observedR) {

		float ret = 0f;
		if (observedR == 0)
			return 0;
		observedR = (float) (observedR * 0.01);
		if (degree == LOH) {

//			if (midreangeLOH(observedR)) {
//				return ret;
//			}
			if (observedR <= 0.5) {

				ret = ((2 * observedR) - 1) / (observedR - 1);

			} else {

				ret = ((2 * observedR) - 1) / observedR;

			}

		} else if (degree == GAIN) {

//			if (midreangeGAIN(observedR)) {
//				return ret;
//			}

			if (observedR < MINForGAIN)
				observedR = MINForGAIN;
			if (observedR > MAXForGAIN)
				observedR = MAXForGAIN;

			if (observedR <= 0.5) {

				ret = (1 - (2 * observedR)) / (observedR);
				//ret = (float) (((double)1/(double)observedR) -2);
				
			} else {

				ret = ((2 * observedR) - 1) / (1 - observedR);

			}

		}
		// System.out.println("tumorratio="+degree+"\t"+observedR+"\t"+ret);
		if((ret > 1) || (ret < 0)){
			ret = 1f;
		}
		return ret;

	}

	private static boolean midreangeLOH(float observedR) {
		// TODO Auto-generated method stub
		return observedR > 0.4 && observedR < 0.6;
	}

	private static boolean midreangeGAIN(float observedR) {
		// TODO Auto-generated method stub
		return observedR > 0.45 && observedR < 0.55;
	}

	public static float getTumorRatioSomatic(float observedS) {
		if (observedS == 0)
			return 0;
		observedS = (float) (observedS * 0.01);
		return observedS*2;
		
	}

}
