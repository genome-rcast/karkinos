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
package jp.ac.utokyo.rcast.karkinos.summary.cnvbygene;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class CNVUtils {

	public static Map<String,GeneBean> getMap(String generef) throws IOException{
		
		//
		Map<String,GeneBean> map = new LinkedHashMap<String,GeneBean>();
		//
		String s = null;
		BufferedReader br = new BufferedReader(new FileReader(generef));
		boolean isFirst = true;
		Map<String,Integer> titles = new HashMap<String,Integer>();
		while ((s = br.readLine()) != null) {
			
			//
			if(isFirst){
				
				isFirst = false;
				String[] t = s.split(",");
				int n = 0;
				for(String ss:t){
					titles.put(ss,n);
					n++;
				}
				continue;
				
			}
			//
			String[] data = s.split(",");
			String geneName = data[titles.get("geneName")];	
			int genelength = Integer.parseInt(data[titles.get("genelength")]);	
			String chrom = data[titles.get("chrom")];	
			int start = Integer.parseInt(data[titles.get("txStart")]);
			int end  = Integer.parseInt(data[titles.get("txEnd")]);
			//
			GeneBean bean = new GeneBean();
			bean.setName(geneName);
			bean.setChr(chrom);
			bean.setStart(start);
			bean.setEnd(end);		
			bean.setLength(genelength);
			
			if(map.containsKey(geneName)){
				
				GeneBean bean0 = map.get(geneName);
				if(bean.getLength()>bean.getLength()){
					map.put(geneName, bean);
				}
				
			}else{
				
				map.put(geneName, bean);
				
			}
			
			
			
		}		
		return map;
		
		
	}
	
	
}
