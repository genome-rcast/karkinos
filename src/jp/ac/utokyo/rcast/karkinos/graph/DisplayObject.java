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
package jp.ac.utokyo.rcast.karkinos.graph;

public class DisplayObject {
	
	public DisplayObject(){
		
	}
	public DisplayObject(Object _obj,int _size,String _title){
		object = _obj;
		size = _size;
		title = _title;
	}
	
	Object object;
	int size;
	public Object getObject() {
		return object;
	}
	public void setObject(Object object) {
		this.object = object;
	}
	public int getSize() {
		return size;
	}
	public void setSize(int size) {
		this.size = size;
	}
	String title ="";
	public String getTitle() {
		// TODO Auto-generated method stub
		return title;
	}
	public void setTitle(String _title) {
		title = _title;		
	}
	
}
