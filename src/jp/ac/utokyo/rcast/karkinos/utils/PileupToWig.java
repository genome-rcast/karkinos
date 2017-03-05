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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class PileupToWig {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String dir ="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/VarScan";
		File f = new File(dir);
		for(File ff:f.listFiles()){
			
			if(ff.getName().endsWith("mpileup")){
				
				String fno = ff.getName();
				String out = fno.replaceAll("mpileup", "wig");
				System.out.println(out);
				outQ20D20(ff,new File(ff.getParentFile().getAbsolutePath()+"/"+out));
				
			}
			
			
		}
		
		

	}

	private static void outQ20D20(File ff, File out) {
		
		//
		try {
		     
			  BufferedWriter bw 
		        = new BufferedWriter(new FileWriter(out, false)); 
		      		      
			
			  int thres = 20;
		      BufferedReader br = new BufferedReader(new FileReader(ff));
		      String line = "";
		      String chrb4 = ""; 
		      int covbf = -1;
		      while ((line = br.readLine()) != null) {

		    	  String[] sa = line.split("\t");
		    	  String chr = sa[0];
		    	  String pos = sa[1];
		    	  int depthn = Integer.parseInt(sa[3]);
		    	  int deptht = Integer.parseInt(sa[6]);
		    	  String qnormal =sa[5];
		    	  String qtumor =sa[8];
		    	  
		    	  int cov = 0;
		    	  if(depthn<20 || deptht<20){
		    		  cov = 0;
		    	  }else{
		    		  cov =1;
		    		  if(depthn<50){
		    			  depthn = depcheck(qnormal);
		    		  }
		    		  if(deptht<50){
		    			  deptht = depcheck(qtumor);
		    		  }
		    		  if(depthn<20 || deptht<20){
			    		  cov = 0;
		    		  }	  
		    	  }		    	  
		    	 if(!chrb4.equals(chr)){
		    		 bw.write("variableStep chrom="+chr+"\n");
		    		 //System.out.println("variableStep chrom="+chr);
		    		 covbf = -1;
		    	 }		    	
		    	//if(covbf != cov){
		    		 bw.write(pos+"\t"+cov+"\n");
		    		 //System.out.println(pos+"\t"+cov);
		    	 //}		    	 
		    	 chrb4 = chr;
		    	 covbf = cov;
		    	 
			    		       
		      }
		      br.close();
		      bw.close();
		    } catch (FileNotFoundException e) {
		     
		      e.printStackTrace();
		    } catch (IOException e) {
		      
		      e.printStackTrace();
		    }
		
		
	}

	private static int depcheck(String s) {
		
		int cnt = 0;
		for(char c:s.toCharArray()){
			
			int phred = c-33;
			
			if(phred>=20){
				cnt++;
			}
		}	
		return cnt;
	}

}
