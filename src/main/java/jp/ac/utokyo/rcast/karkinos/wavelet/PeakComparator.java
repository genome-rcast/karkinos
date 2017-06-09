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

import java.util.Comparator;
import java.util.List;

public class PeakComparator implements Comparator<Peak> {

	public int compare(Peak p1, Peak p2) {
		
		if(p1.getU()==p2.getU())return 0;
		return p1.getU()<p2.getU()?-1:1;
	}

}
