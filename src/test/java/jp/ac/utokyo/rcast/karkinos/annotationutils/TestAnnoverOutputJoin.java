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
package jp.ac.utokyo.rcast.karkinos.annotationutils;

import java.util.ArrayList;
import java.util.List;

public class TestAnnoverOutputJoin {
	
	
	public static void main(String[] arg){
		
		//
		String s1="/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_baylor/HCC-JP-75-T-HCC-JP-75-N/HCC-JP-75-T-HCC-JP-75-N" +
				"/HCC-JP-75-T-HCC-JP-75-N_annover_input.txt.genome_summary.csv";
		String s2="/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_baylor/HCC-JP-75-T-HCC-JP-75-N/HCC-JP-75-T-HCC-JP-75-N/HCC-JP-75-T-HCC-JP-75-N_snvdata.vcf";
		List<String> l = new ArrayList<String>();		
		add(l,"-a",s1);
		add(l,"-v",s2);
		String[] ar = l.toArray(new String[l.size()]);
		AnnoverOutputJoin.main(ar);
		
	}
	
	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}
}
