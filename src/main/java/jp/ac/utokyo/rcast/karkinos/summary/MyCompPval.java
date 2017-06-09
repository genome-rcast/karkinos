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
package jp.ac.utokyo.rcast.karkinos.summary;

import java.util.Comparator;
import java.util.List;

public class MyCompPval implements Comparator<List> {

	int idx = 0;
	public MyCompPval(int i) {
		idx = i;
	}

	public int compare(List arg0, List arg1) {
		
		Object pval0 = arg0.get(idx);
		Object pval1 = arg1.get(idx);
		double d0 = Double.parseDouble(pval0.toString());
		double d1 = Double.parseDouble(pval1.toString());
		if(d0!=d1){
		 return (d0>d1)?1:-1;
		}
		Object lval0 = arg0.get(idx+1);
		Object lval1 = arg1.get(idx+1);
		double l0 = Double.parseDouble(lval0.toString());
		double l1 = Double.parseDouble(lval1.toString());
		if(l0==l1)return 0;
		return (l0>l1)?-1:1;
	}
	
}
