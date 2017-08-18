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
package jp.ac.utokyo.rcast.karkinos.annotation;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class TestRead {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		// 
		String cosmic  = "/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf";
		FileInputStream fis = null;
		BufferedReader br = null;
		try {
			fis = new FileInputStream(cosmic);			   
			
			fis.skip(90896);
			br = new BufferedReader(new InputStreamReader(fis));
			int totalcnt = 0;
			long init = fis.available();
			int loadtotal = 0;
			int cntignore = 0;
			String line = br.readLine();
			System.out.println(line);	
			line = br.readLine();
			System.out.println(line);	
		} catch (Exception ex) {
				
		}		
		
	}
	
     public static void _main(String[] args) {
		
		// 
		String cosmic  = "/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf";
		FileInputStream fis = null;
		BufferedReader br = null;
		try {
			fis = new FileInputStream(cosmic);			   
			
			fis.skip(30711);
			InputStreamReader is = (new InputStreamReader(fis));
			int totalcnt = 0;
			long init = fis.available();
			int loadtotal = 0;
			int cntignore = 0;
			String line = br.readLine();
			System.out.println(line);	
			line = br.readLine();
			System.out.println(line);	
		} catch (Exception ex) {
				
		}		
		
	}

}
