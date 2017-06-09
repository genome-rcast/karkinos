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
package jp.ac.utokyo.rcast.karkinos.annotation.loadsave;

import java.util.ArrayList;
import java.util.List;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;

public class LoadUtils {

	public static void merge(SaveBean sb, SaveBean sbToAdd){
		
		//
		DataSet ds1 = sb.getDataset();
		DataSet ds2 = sbToAdd.getDataset();
		
		ReadsSummary rs1 = sb.getReadsSummary();
		ReadsSummary rs2 = sbToAdd.getReadsSummary();
		
		DataSet dataset = merge(ds1,ds2);
		ReadsSummary readsSummary = merge(rs1,rs2);
		
		sb.setDataset(dataset);
		sb.setReadsSummary(readsSummary);
		
	}



	private static DataSet merge(DataSet ds1, DataSet ds2) {
		
		if(ds1==null)return ds2;
		addall(ds1.getSnvlist(),ds2.getSnvlist());
		ds1.getNormal().merge(ds2.getNormal());
		ds1.getTumor().merge(ds2.getTumor());
		
		return ds1;
	}
	
	private static void addall(List list, List list2) {
		if(list!=null&&list2!=null){
			list.addAll(list2);
		}else if(list!=null){
			//donothing
		}else{
			list = new ArrayList();
			list.addAll(list2);
		}
	}



	private static ReadsSummary merge(ReadsSummary rs1, ReadsSummary rs2) {
		
		if(rs1==null)return rs2;
		//
		rs1.merge(rs2);		
		return rs1;
	}

}
