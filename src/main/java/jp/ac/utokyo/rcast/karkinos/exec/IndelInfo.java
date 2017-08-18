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

public class IndelInfo implements java.io.Serializable{
	
	public boolean indel = false;
	public String insersion;
	public int length;
	public int cnt;
	public boolean reg = false;
	public int refpos;;
	public void clear() {
		reg = true;
		indel = false;
		insersion=null;
		length = 0;
	}
	public boolean isIndel() {
		return indel;
	}
	public String getInsersion() {
		return insersion;
	}
	public int getLength() {
		return length;
	}
	public int getCnt() {
		return cnt;
	}	
	
}
