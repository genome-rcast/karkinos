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
package jp.ac.utokyo.rcast.karkinos.annotation.stats;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

public class Fisher {
	
	
	public static void main(String[] arg){
		
		System.out.println(calcPValue(8,2,8,0));
		
		
	}
	
	
	public static double calcPValue(int a,int b,int c,int d) {
		
		int[][] data = new int[2][2];
		data[0][0] = a;
		data[0][1] = b;
		data[1][0] = c;
		data[1][1] = d;
		return calcPValue(data);
	}
	
	/**
	 * Pvalue fisher
	 *
	 * @param data
	 *            ２×２matrix
	 * @return P値
	 */
	public static double calcPValue(int[][] data) {
		checkInputData(data);

		double pValue = 0;
		double criterion = calcProbability(data);
		List feasibleTableList = constractFeasibleTableList(data);

		for (Iterator iter = feasibleTableList.iterator(); iter.hasNext();) {
			int[][] dataTmp = (int[][]) iter.next();
			double p = calcProbability(dataTmp);

			if (p <= criterion) {
				pValue += p;
			}
		}

		if (pValue > 1) {
			return 1;
		} else {
			return pValue;
		}
	}

	
	
	/**
	 * feasible list
	 */
	private static List constractFeasibleTableList(int[][] data) {
		checkInputData(data);

		int a = data[0][0];
		int b = data[0][1];
		int c = data[1][0];
		int d = data[1][1];

		int e = a + b;
		int g = a + c;
		int h = b + d;

		List feasibleTableList = new ArrayList();
		for (int a2 = 0; a2 <= Math.min(e, g); a2++) {
			int b2 = e - a2;
			int c2 = g - a2;
			int d2 = h - b2;

			if (b2 < 0 || c2 < 0 || d2 < 0) {
				continue;
			}

			feasibleTableList.add(new int[][] { { a2, b2 }, { c2, d2 } });
		}
		return feasibleTableList;
	}
	/**
	 * result = ((a+b)!/a! (b+d)!/b! (a+c)!/c! (c+d)!/d!) / (a+b+c+d)!
	 */
	private static double calcProbability(int[][] data) {
		checkInputData(data);

		int a = data[0][0];
		int b = data[0][1];
		int c = data[1][0];
		int d = data[1][1];

		int n = 1;
		int max = Math.max(a + b, Math.max(b + d, Math.max(a + c, c + d)));

		double result = 1;
		/*calc from small val */
		for (int i = 1; i <= max; i++) {
			if (a < i && i <= a + b) { /* (a+b)! / a! */
				result *= (double) i / n;
				n++;
			}
			if (b < i && i <= b + d) { /* (b+d)! / b! */
				result *= (double) i / n;
				n++;
			}
			if (c < i && i <= a + c) { /* (a+c)! / c! */
				result *= (double) i / n;
				n++;
			}
			if (d < i && i <= c + d) { /* (c+d)! / d! */
				result *= (double) i / n;
				n++;
			}
		}
		return result;
	}
	/**
	 * 
	 */
	private static void checkInputData(int[][] data) {
		if (data.length != 2 || data[0].length != 2) {
			throw new IllegalArgumentException();
		}

		int a = data[0][0];
		int b = data[0][1];
		int c = data[1][0];
		int d = data[1][1];

		if (a < 0 || b < 0 || c < 0 || d < 0) {
			throw new IllegalArgumentException();
		}
	}

	private static int max(int n1,int n2){
		return Math.max(n1, n2);
	}
	
	/**
	 * オッズ比を求める
	 *
	 * @param data
	 *            ２×２の非負整数行列
	 * @return オッズ比
	 */
	public static double calcOddsRatio(int _a, int _b,int _c,int _d) {
		
		double a = (double) _a;
		double b = (double) _b;
		double c = (double) _c;
		double d = (double) _d;


		if (a == 0 || b == 0 || c == 0 || d == 0) {
			a += 0.5;
			b += 0.5;
			c += 0.5;
			d += 0.5;
		}
		double result = (a * d) / (b * c);

		return result;
	}

	/**
	 * オッズ比の下側95%信頼限界値を求めます。
	 *
	 * @param data
	 *            ２×２の非負整数行列
	 * @return オッズ比の下側95%信頼限界値
	 */
	public static double calcOddsRatioCI(int _a, int _b,int _c,int _d,boolean posi,double bound) {

		double a = (double) _a;
		double b = (double) _b;
		double c = (double) _c;
		double d = (double) _d;

		if (a == 0 || b == 0 || c == 0 || d == 0) {
			a += 0.5;
			b += 0.5;
			c += 0.5;
			d += 0.5;
		}
		double oddsRatio = calcOddsRatio(_a,_b,_c,_d);
		if(!posi){
			bound = bound * -1;
		} 
		double result = oddsRatio
				* Math.exp(bound
						* Math.sqrt(1 / a + 1 / b + 1 / c + 1 / d));
		return result;
	}


}
