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

import java.awt.geom.Point2D;
import java.util.List;

public class RegressionInfo implements java.io.Serializable {

	List<Point2D> gcMedianList = null;
	// y=ae^bx+c
	double a = 0, b = 0, c = 0, d=0;
	double baseAdjust = 0;
	public double getBaseAdjust() {
		return baseAdjust;
	}

	public void setBaseAdjust(double baseAdjust) {
		this.baseAdjust = baseAdjust;
	}

	public double getD() {
		return d;
	}

	public void setD(double d) {
		this.d = d;
	}

	int funcflg = 0;


	public int getFuncflg() {
		return funcflg;
	}

	public void setFuncflg(int funcflg) {
		this.funcflg = funcflg;
	}

	public List<Point2D> getGcMedianList() {
		return gcMedianList;
	}

	public void setGcMedianList(List<Point2D> gcMedianList) {
		this.gcMedianList = gcMedianList;
	}

	public double getA() {
		return a;
	}

	public void setA(double a) {
		this.a = a;
	}

	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}

	public double getC() {
		return c;
	}

	public void setC(double c) {
		this.c = c;
	}

	public double getReg(double x) {
		
		if (funcflg==FunctionRegression.Exp) {
			double y = a * Math.exp(b * x) + c;
			return y+baseAdjust;
		} else if(funcflg==FunctionRegression.Log) {
			double y = a * Math.log(x) + b;
			return y+baseAdjust;
		}else if(funcflg==FunctionRegression.Poly2){
			double y = a * Math.pow(x, 2)+b*x+c;
			return y+baseAdjust;
		}else if(funcflg==FunctionRegression.Poly3){
			double y = a * Math.pow(x, 3)+b * Math.pow(x, 2)+c*x+d;
			return y+baseAdjust;
		}
		return 0;
	}

	public double getRegOrg(double x) {
		if (funcflg==FunctionRegression.Exp) {
			double y = a * Math.exp(b * x) + c;
			return y;
		} else if(funcflg==FunctionRegression.Log) {
			double y = a * Math.log(x) + b;
			return y;
		}else if(funcflg==FunctionRegression.Poly2){
			double y = a * Math.pow(x, 2)+b*x+c;
			return y;
		}else if(funcflg==FunctionRegression.Poly3){
			double y = a * Math.pow(x, 3)+b * Math.pow(x, 2)+c*x+d;
			return y;
		}
		return 0;
	}

}
