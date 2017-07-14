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
package jp.ac.utokyo.rcast.karkinos.cntwavelet;

import java.util.ArrayList;
import java.util.List;

public class MexicanHatWavelet {

	public static void main(String[] arg) {

		//
		double data[] = new double[4000];
		int idx  =0;
		for(double d:data){
			data[idx]= 1;
			idx++;
		}
		CTWaveletBean bean = getWaveletTransform(data);
		
//		double sum = 0;
//		for(double d= -1;d<1;d=d+0.001){
//			
//			double v = getMexicanHat(0.1,d);
//			System.out.println(d +"\t"+ v);
//			sum = sum+v;
//		}
//		System.out.println(sum);
	}
	
	public static CTWaveletBean getWaveletTransform(double[] data) {

		//
		CTWaveletBean bean = new CTWaveletBean();
		tf(bean, data);
		return bean;

	}

	private static void tf(CTWaveletBean bean, double[] data) {

		//
		// start sd = 0.15 to 0.01
		double sumd = 0;
		for(double dd:data){
			sumd = sumd+dd;
		}
		int l=0;
		for(double d:data){
			data[l] = d/sumd;
			l++;
		}
		
		double sds = 0.01;
		double sd = sds;
		double sde = 0.03;
		double sdintv = 0.001;

		double maxpeak = -10;
		double maxsd = 0;
		double[] maxfit = null;
		double[] data2 = null;

		
		int m = 0;
		
		for (sd = sds; sd < sde; sd = sd + sdintv) {

			//			
			data2 = new double[data.length];	
			double sumval = 0;
			for (int n = 0; n < 3999; n = n + 1) {
				//
				double val = tf(data, sd, n);
				//System.out.println(val);
				data2[n] = val;
				sumval   =sumval +val;
//				if(val  > maxpeak){
//					maxpeak = val ;
//					maxfit = data2;
//					maxsd = sd;
//					double mostfittingvariance = Math.pow(sd, 2);
//					bean.setMostfittingvariance(mostfittingvariance);
//					System.out.println("val="+val+"\t"+sd);
//				}

			}

			if (sumval  > maxpeak) {
				maxpeak = sumval ;
				maxfit = data2;
				maxsd = sd;
				double mostfittingvariance = Math.pow(sd, 2);
				bean.setMostfittingvariance(mostfittingvariance);
				System.out.println("sum="+sumval+"\t"+sd);
				bean.setData(data2);
				bean.setSd(sd);
				sumval = 0;
			}
			
			m++;

		}
		
		

	}

	private static double tf(double[] data, double sd, int n) {

		int start = n - 300;
		int end = n + 300;

		double x = (n * 0.001);
		double xs = (start * 0.001);
		if (xs < 0) {
			xs = 0;
			start = 0;
		}

		double xe = (end * 0.001);
		if (xe > 4) {
			xe = 4;
			end = 3999;
		}
		//
		double sum = 0;
		for (double d = xs; d <= xe; d = d + 0.001) {
			//
			int idx = (int) (d * 1000);
			if (idx > data.length - 1) {
				idx = data.length - 1;
			}
			//
			if(start>=data.length){
				break;
			}
			double val = data[idx];
			double diff = d - x;
			double v = getMexicanHat(sd, diff);
			
			double product = (v * val);
			//System.out.println("product="+product+"\t"+idx+"\t"+x);
			sum = sum + product;

		}
		return sum;
	}

	private static double getMexicanHat(double sd, double diff) {

		
		double x = (diff / sd);
		double x2 = pow2(x);
		double c = 2.0 / cube(Math.PI);
		double f1 = (1.0 / sqrt(3 * (sd)));
		double f2 = (1 - x2);
		double ev = Math.exp(-0.5 * x2);
		return (c* f1 * f2 * ev);

	}
	
//	private static double getMexicanHat(double sd, double diff) {
//
//		diff = diff * sd;
//		
//		double x = (diff);
//		double x2 = pow2(x);
//		double c = 2.0 / cube(Math.PI);
//		double f1 = (1.0 / sqrt(3));
//		double f2 = (1 - x2);
//		double ev = Math.exp(-0.5 * x2);
//		return (c* f1 * f2 * ev);
//
//	}

	public static double pow2(double d) {
		return Math.pow(d, 2);
	}

	public static double sqrt(double d) {
		return Math.sqrt(d);
	}

	public static double cube(double d) {
		return sqrt(sqrt(d));
	}

}
