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
package jp.ac.utokyo.karkinos.noisefilter;

public class NDist {
	
	double u  =0;
	double v =0;
	
	public NDist(double u,double v){
		
		//
		this.u =u;
		this.v =v;		
	}
	
	public double getP(double x){
		
		double p = Math.exp(-0.5*(Math.pow((x-u),2)/v))/Math.sqrt(2*Math.PI*v);
		return p;
	}

}
