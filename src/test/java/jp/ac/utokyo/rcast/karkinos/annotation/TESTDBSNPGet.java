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
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CaptureHolder;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

public class TESTDBSNPGet {

	public static void main(String[] arg){
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_snp132.txt";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_ALL.sites.2012_02.txt";
		String cosmic  = "/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf";
		DbSNPAnnotation dbAnno = new DbSNPAnnotation(dbSNP, g1000,0.01f,cosmic,null);
		DbSNPBean snpBean =  dbAnno.lookup("1", 866408);
		System.out.println(snpBean.getInfo());
		
	}
	
	public static void _main(String[] arg) throws IOException{
		
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_snp132.txt";
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_ALL.sites.2012_02.txt";
		
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		
		CaptureHolder ch = new CaptureHolder();
		ch.loadTargetBed(targetRegion,tgr);
		//
		
		FileInputStream fis = new FileInputStream(g1000);
		BufferedReader br = new BufferedReader(new InputStreamReader(fis));

		int cntsnp = 0;
		int cntsnpgene = 0;
		int cntsnpnc = 0;
		
		try {

			
			int totalcnt = 0;
			int init = fis.available();
			String chr = "";
			for (;;) {
				int point = init - fis.available();
				String line = br.readLine();
				if (line == null)
					break;
				totalcnt++;
				String[] sa = line.split("\t");
				String _chr = sa[1];
				if (!chr.equals(_chr)) {
					System.out.println(_chr + "\t" + totalcnt + "\t" + point);
				}
				chr = _chr;
				int pos = Integer.parseInt(sa[2]);
				
				CapInterval ci = ch.getCapInterval(_chr, pos);
				if(ci!=null){
					cntsnp++;
					if(ci.isGene()){
						cntsnpgene++;
					}else{
						cntsnpnc++;
					}
					
				}
				
			}

		} finally {
			br.close();
		}		
		System.out.println("cnt snp ="+cntsnp);
		System.out.println("cnt snp gene ="+cntsnpgene);
		System.out.println("cnt snp nc ="+cntsnpnc);
	}
	
}
