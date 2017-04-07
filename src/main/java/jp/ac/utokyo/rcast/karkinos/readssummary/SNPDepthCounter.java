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
package jp.ac.utokyo.rcast.karkinos.readssummary;

import jp.ac.utokyo.rcast.karkinos.exec.KarkinosProp;

public class SNPDepthCounter {
	
	
	int lowdepthcnt=0;
	int x2count=0;
	int x3count=0;
	int x4count=0;
	int highdepthcnt=0;
	int total = 0;
	
	public void reg(int depth){
		
		if(depth > (KarkinosProp.mindepth*5)){
			highdepthcnt++;	
		}else if(depth > (KarkinosProp.mindepth*4)){
			x4count++;	
		}else if(depth > (KarkinosProp.mindepth*3)){
			x3count++;	
		}else if(depth > (KarkinosProp.mindepth*2)){
			x2count++;		
		}else{
			lowdepthcnt++;
		}
		total++;
		
	}

	public int getTotal() {
		return total;
	}

	public int getLowdepthcnt() {
		return lowdepthcnt;
	}

	public int getX2count() {
		return x2count;
	}

	public int getX3count() {
		return x3count;
	}

	public int getX4count() {
		return x4count;
	}

	public int getHighdepthcnt() {
		return highdepthcnt;
	}

}
