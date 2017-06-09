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
package jp.ac.utokyo.rcast.karkinos.utils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class PositionBrocks {
		
	TreeMap<Integer,Integer> treemap 
		= new TreeMap<Integer,Integer>();
	
	public void read(int count,LittleEndian la) throws IOException{
		
		List<Integer> intv = new ArrayList<Integer>();
		for(int n=0;n<count;n++)intv.add(la.readInt());
		List<Integer> sizes = new ArrayList<Integer>();
		for(int n=0;n<count;n++)sizes.add(la.readInt());
		
		int idx = 0;
		for(int pos:intv){
			
			int size = sizes.get(idx);
			treemap.put(pos, size);
			//System.out.println(pos+":"+size);
			idx++;
		}
		intv  = null;
		sizes = null;
		
	}
	
	public boolean inBlock(int pos){
		
		if(treemap==null)return false;
		if(treemap.containsKey(pos)){
			return true;
		}
		Integer floorPos = treemap.floorKey(pos);
		if(floorPos==null)return false;
		int floorBlockSize = treemap.get(floorPos);
		int end = (floorPos+floorBlockSize);
		return (pos<end);
		
	}
	
	
}
