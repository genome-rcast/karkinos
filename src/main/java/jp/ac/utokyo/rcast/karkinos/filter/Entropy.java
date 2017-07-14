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
package jp.ac.utokyo.rcast.karkinos.filter;

public class Entropy {

	public static void main(String[] arg){
		
		System.out.println(entropy("Ttaaaaaaaaaa"));
		System.out.println(entropy("ATGC"));
	}
	public static double entropy(String seq){
		
		
		seq = seq.toUpperCase();
		double[] rcounter = countRatio(seq);
		double s = 0;
		for(double d:rcounter){
			
			if(d==0)continue;
			double sn = -1*d*log2(d);
			s = s+ sn;
		}
		return s;
		
	}

	private static double log2(double d) {
		
		return Math.log(d)/Math.log(2);
	}

	private static double[] countRatio(String s) {
		
		double[] da = new double[4];
		int a=0,t=0,g=0,c=0,n=0;
		int len = s.length();
		for(char ca:s.toCharArray()){
			
			if(ca=='A'){
				a++;
			}else if(ca=='T'){
				t++;
			}else if(ca=='C'){
				c++;
			}else if(ca=='G'){
				g++;
			}else if(ca=='N'){
				n++;
			}
			
		}
		len = len-n;
		da[0]=(double)a/(double)len;
		da[1]=(double)t/(double)len;
		da[2]=(double)c/(double)len;
		da[3]=(double)g/(double)len;
		
		return da;
	}
	
}
