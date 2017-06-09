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

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class GCParcentAdjust {

	public static void calc(DataSet dataset) {

		List<double[]> xyzList = new ArrayList<double[]>();
		List<List<WaveletIF>> cap = dataset.getCapInterval();
		int minbaitlen = 0;
		for (List<WaveletIF> list : cap) {

			for (WaveletIF wi : list) {

				CapInterval civ = (CapInterval) wi;
				double[] ar = new double[3];
				ar[0] = civ.getLength();
				ar[1] = civ.getCgParcent();
				ar[2] = civ.getOriginalValue();
				if(minbaitlen==0||minbaitlen>civ.getLength()){
					minbaitlen = civ.getLength();
				}
				if(civ.getLength()>0){
					xyzList.add(ar);
				}
			}

		}

		//
		FunctionRegression fr = new FunctionRegression(xyzList,minbaitlen);
		dataset.setGcFunctionRegression(fr);
		
		for (List<WaveletIF> list : cap) {

			for (WaveletIF wi : list) {

				CapInterval civ = (CapInterval) wi;
				
				double x = civ.getLength();
				double y = civ.getCgParcent();
				double z = civ.getOriginalValue();
				//
				double adjustedZ = fr.getAdjustedZ(x, y, z);
				civ.setGCAdjustedTNratio(adjustedZ);
				
			}

		}

	}
}
