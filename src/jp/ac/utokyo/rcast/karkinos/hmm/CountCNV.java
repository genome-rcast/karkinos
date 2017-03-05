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
package jp.ac.utokyo.rcast.karkinos.hmm;

import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class CountCNV {

	public static int count(DataSet dataset) {
		int cnt = 0;
		List<List<WaveletIF>> plist = dataset.getCapInterval();
		int b4 = 0;
		for (List<WaveletIF> list : plist) {

			for (WaveletIF wi : list) {
				
				CapInterval ci = (CapInterval)wi;
				int cn = ci.getPeakIdx();
				if(cn!=b4){
					cnt++;
				}				
				b4 = cn;
			}
	
		}
		return cnt/2;
	}

}
