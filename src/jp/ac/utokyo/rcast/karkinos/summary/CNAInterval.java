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
package jp.ac.utokyo.rcast.karkinos.summary;

public class CNAInterval {

	public String getChr() {
		return chr;
	}

	public String getName() {
		return name;
	}

	//bw.write("#chr"+"\t");
	//bw.write("start"+"\t");
	//bw.write("end"+"\t");
	//bw.write("copy number"+"\t");
	//bw.write("gain/loss"+"\t");
	//if(dispFlg==1){
	// bw.write("num hetro SNP"+"\t");
	// bw.write("varidated by allelic"+"\t");
	//}
	public CNAInterval(String[] data) {
		chr = data[0];
		start = Integer.parseInt(data[1].replaceAll(",", ""));
		end = Integer.parseInt(data[2].replaceAll(",", ""));
		try{
		 cn = Float.parseFloat(data[3].replaceAll("n=", ""));
		}catch(Exception ex){
			System.out.println(data);
		}
		if (data.length > 7) {
			try {
				noSNP = Integer.parseInt(data[5].replaceAll(",", ""));
				//snpclrrel = Float.parseFloat(data[6]);
				aaf = Float.parseFloat(data[6]);
				baf = Float.parseFloat(data[7]);
				//varidated = data[6];
				type = data[8];
			} catch (Exception ex) {
			}
		}else{
			type = data[5];
		}
	}

	public CNAInterval() {

	}

	String chr;
	int start;
	int end;
	float cn;
	String name;
	private String gieStain;
	String varidated;
	float snpclrrel;
	int noSNP;
	String type;
	float aaf;
	float baf;
	public boolean intercect(int start2, int end2) {
		
		return start<end2 && start2 < end;
	}

	public void setGieStain(String gieStain) {
		this.gieStain = gieStain;
	}

	public String getGieStain() {
		return gieStain;
	}
	
}
