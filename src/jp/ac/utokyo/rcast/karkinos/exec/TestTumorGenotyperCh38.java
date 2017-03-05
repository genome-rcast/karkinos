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

import java.util.ArrayList;
import java.util.List;

public class TestTumorGenotyperCh38 {
	
	

	
	public static void main(String[] arg){
		
		
		String normalbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3N_TDv3_genome.bam";
		String tumorbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3TA_TDv3_genome.bam";
		
		String twobitref = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/hg38.2bit";
		//String dbSNP = "/GLUSTER_DIST/data/users/ueda/SNVtest/hg19_ALL.sites.2012_02.txt";
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/snp147Common_hg38_ontarget.txt";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/1000g_phase3_ontarget_trim.vcf";
		
		String cosmic = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/CosmicCodingMuts.vcf";
		String refflat = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/refFlat_panel_sort.txt";
		
		//String exomesnp = "/usr/local/karkinos/karkinos/genome/exomeSNP.vcf";
		
		//String targetRegion = "/data/users/ueda/project/TodaiPanel/S3035822_Covered.bed_test";
		String targetRegion = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/S3035822_Covered_sm_lift_hg38.bed";
		
		//String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/vcrome2.1.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/halo_ver2_Regions.bed";
		//String mappability = "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";
		//String mappability = "/GLUSTER_DIST/data/Genomes/karkinos/genome/wgEncodeCrgMapabilityAlign36mer.bw";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinostest2/";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/project/gotoh/S1-S7";
		
		//String outdir = "/GLUSTER_DIST/data/users/ueda/test2";
		String outdir = "/GLUSTER_DIST/data/users/ueda/toptesthg38";
		String propfile = "/usr/local/karkinos/karkinos4.0/karkinos.properties";
		
		String definedList = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/hotspot.txt";
		
		//
		List<String> l = new ArrayList<String>();
		//add(l,"-n",normalbamf);
		//add(l,"-md",middelfile);
		add(l,"-n",normalbamf);
		add(l,"-t",tumorbamf);
		//add(l,"-chr","chr1");
		//add(l,"-startend","2");
		add(l,"-r",twobitref);
		add(l,"-snp",dbSNP);
		add(l,"-ct",targetRegion);
		//add(l,"-mp",mappability);
		//add(l,"-id","CRC_107_T-CRC_107_N");
		add(l,"-id","TK-3");
		add(l,"-o",outdir);
		add(l,"-prop",propfile);
		add(l,"-g1000",g1000);
		add(l,"-g1000freq","0.01");
		//add(l,"-tc","1.0");
		//add(l,"-rg","/GLUSTER_DIST/data/Genomes/ucscgene_hg19/totoki_hg19/refFlat_120402.txt");
		
		add(l,"-rg",refflat);
		add(l,"-cosmic",cosmic);
		//add(l,"-exonSNP",exomesnp);
		
		add(l,"-sites",definedList); 
				
		//add(l,"-nd","true");
		
		//add(l,"-rs","/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/total_reads_stats_n.txt," +
		//		"/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/total_reads_stats_t.txt");
		
		String[] ar = l.toArray(new String[l.size()]);
		//TumorGenotyperReanalysis.main(ar);
		TumorGenotyper.main(ar);
		
	}

	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

}
