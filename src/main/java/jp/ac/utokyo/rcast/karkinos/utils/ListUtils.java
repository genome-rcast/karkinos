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

import java.util.ArrayList;
import java.util.List;

public class ListUtils {

	public static List<Interval> getIntervalList(String chrom, String startend,
			int length, int binbitsize) {

		//
		int s = 0;
		int e = 0;
		if (startend != null) {
			String[] sa = startend.split("-");
			
			if(sa.length>1){
				s = Integer.parseInt(sa[0]);
				e = Integer.parseInt(sa[1]);
			}else{
				int chank = Integer.parseInt(sa[0]);
				int unit = 20000000;
				s = (chank-1)*(unit);
				e = (chank)*(unit);
				if(s>length){
					return null;
				}
			}
			
		}

		List<Interval> list = new ArrayList<Interval>();
		int binsize = (int) Math.pow(2, binbitsize);
		int start = 0;
		int end = binsize;
		if (s != 0) {
			start = s;
			end = start+binsize;
		}
		if (e != 0 && e < length) {
			length = e;
		}	
		while (end < length) {

			Interval iv = new Interval(chrom, start, end);
			list.add(iv);
			start = end + 1;
			end = end + binsize;
			
		}
		Interval iv = new Interval(chrom, start, length);
		list.add(iv);
		return list;

	}

}
