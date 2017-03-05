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

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;

public class SaveBean implements java.io.Serializable{
	

	private static final long serialVersionUID = -3177838380080434442L;


	public SaveBean(DataSet _dataset,ReadsSummary _readsSummary){
			
		//
		dataset = _dataset;
		readsSummary = _readsSummary;
	}
	
	ReadsSummary readsSummary;
	DataSet dataset;
	
	public ReadsSummary getReadsSummary() {
		return readsSummary;
	}
	public void setReadsSummary(ReadsSummary readsSummary) {
		this.readsSummary = readsSummary;
	}
	public DataSet getDataset() {
		return dataset;
	}
	public void setDataset(DataSet dataset) {
		this.dataset = dataset;
	}
	
}
