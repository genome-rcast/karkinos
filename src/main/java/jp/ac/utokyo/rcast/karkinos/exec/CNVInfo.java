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
package jp.ac.utokyo.rcast.karkinos.exec;

public class CNVInfo implements java.io.Serializable{
	
	long normalcnt;
	double normaldepth;
	double normaldepthAdj;
	long tumorcnt;
	double tumordepth;
	double tumordepthAdj;
	double tnratio;
	double denoise;
	double copynumber;
	double tnratioorg;
	
	public CNVInfo(long _normalcnt, double _normaldepth,
			double _normaldepthAdj, long _tumorcnt, double _tumordepth,
			double _tumordepthAdj, double _tnratio) {
		
		normalcnt = _normalcnt;
		normaldepth = _normaldepth;
		normaldepthAdj = _normaldepthAdj;
		tumorcnt = _tumorcnt;
		tumordepth = _tumordepth;
		tumordepthAdj = _tumordepthAdj;
		tnratio =_tnratio;
		tnratioorg = _tnratio;
	}
	public long getNormalcnt() {
		return normalcnt;
	}

	public double getNormaldepth() {
		return normaldepth;
	}

	public double getNormaldepthAdj() {
		return normaldepthAdj;
	}

	public long getTumorcnt() {
		return tumorcnt;
	}

	public double getTumordepth() {
		return tumordepth;
	}

	public double getTumordepthAdj() {
		return tumordepthAdj;
	}

	public double getTnratio() {
		return tnratio;
	}
	public void setDenioseValue(double _denoise) {
		denoise = _denoise;
	}
	public void setCN(double _copynumber) {
		copynumber = _copynumber;
	}
	public double getDenoise() {
		return denoise;
	}
	public double getCopynumber() {
		return copynumber;
	}
	public double getOriginalTnratio() {
		return tnratioorg;
	}
	

}
