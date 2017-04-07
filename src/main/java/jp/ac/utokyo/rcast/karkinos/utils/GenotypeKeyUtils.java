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

public class GenotypeKeyUtils {
	
	public static final String keys1[] = new String[] { "CtoT", "CtoG", "CtoA",
		"TtoG", "TtoC", "TtoA" };
	public static final String keys2[] = new String[] { "GtoA", "GtoC", "GtoT",
		"AtoC", "AtoG", "AtoT" };
	
	public static String aggrigateKeys(String key) {
		if (index(keys1, key) >= 0) {
			return key;
		} else {
			int idx = index(keys2, key);
			if (idx == -1) {
				return "";
			}
			return keys1[idx];
		}

	}	
	
	public static String toDispKey(String key){
		int idx1  = index(keys1, key);
		String key2 = keys2[idx1];
		return key.charAt(0)+":"+key2.charAt(0)+">"
				+key.charAt(3)+":"+key2.charAt(3);
	}
	
	public static int index(String[] keys, String key) {

		int n = 0;
		for (String s : keys) {
			if (s.equals(key)) {
				return n;
			}
			n++;
		}
		return -1;
	}

	
}
