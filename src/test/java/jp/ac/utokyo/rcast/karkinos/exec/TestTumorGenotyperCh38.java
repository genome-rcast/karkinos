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
		
		
		
		//String normalbamf ="/home/spark/todaitoptest/testdata0308/similate/LUAD-311-N_TDv3_genome_0.bam";
		//String tumorbamf ="/home/spark/todaitoptest/testdata0308/similate/LUAD-311-N_TDv3_genome_1_2.bam";

		
		//String normalbamf = "/home/spark/todaitoptest/testdata0308/realign/LUAD-311-N_TDv3_genome.bam_realign.bam";
		//String tumorbamf = "/home/spark/todaitoptest/testdata0308/realign/LUAD-311-T_TDv3_genome.bam_realign.bam";
		
		System.out.println("test start");
		
		//String normalbamf = "/home/spark/todaitoptest/test/664-9136-6-N_TDv3-NS_H17-3266-T_TDv3-NS/realign/664-9136-6-N_TDv3-NS_bwa.sorted.dedup.bam_realign.bam";
		//String tumorbamf = "/home/spark/todaitoptest/test/664-9136-6-N_TDv3-NS_H17-3266-T_TDv3-NS/realign/H17-3266-T_TDv3-NS_bwa.sorted.dedup.bam_realign.bam";
		

		String normalbamf = "/home/spark/todaitoptest/testdata/7564429-N_TDv3-NS_bwa.sorted.dedup.bam_realign.bam";
		String tumorbamf = "/home/spark/todaitoptest/testdata/H15-14183-23-T_TDv3-NS_bwa.sorted.dedup.bam_realign.bam";
		
		//String normalbamf ="/home/spark/todaitoptest/testdata0308/similate/LUAD-311-N_TDv3_genome_1.bam";
		//String tumorbamf ="";
				
		//String normalbamf = "/home/spark/todaitoptest/testdata0308/LUAD-311-N_TDv3_genome.bam";
		//String tumorbamf = "/home/spark/todaitoptest/testdata0308/LUAD-311-T_TDv3_genome.bam";
		
//		String normalbamf ="/home/spark/todaitoptest/testdata0308/LUAD-311-N_TDv3_genome.bam";
//		String tumorbamf = "/home/spark/todaitoptest/testdata0308/LUAD-311-T_TDv3_genome.bam";
		
//		String normalbamf = "/home/spark/todaitoptest/testdata/LUAD-311N_0_2_4sort.bam";
//		String tumorbamf = "/home/spark/todaitoptest/testdata/LUAD-311N_6_12sort.bam";
		
		String twobitref = "/home/spark/todaitoptest/ref/hg38.2bit";
		//String dbSNP = "/GLUSTER_DIST/data/users/ueda/SNVtest/hg19_ALL.sites.2012_02.txt";
		String dbSNP = "/home/spark/todaitoptest/ref/snp147Common_hg38_ontarget.txt";
		String g1000 = "/home/spark/todaitoptest/ref/1000g_phase3_ontarget_trim.vcf";
		
		String cosmic = "/home/spark/todaitoptest/ref/CosmicCodingMuts.vcf";
		String refflat = "/home/spark/todaitoptest/ref/refFlat_panel_sort.txt";
		
		//String exomesnp = "/usr/local/karkinos/karkinos/genome/exomeSNP.vcf";
		
		//String targetRegion = "/data/users/ueda/project/TodaiPanel/S3035822_Covered.bed_test";
		String targetRegion = "/home/spark/todaitoptest/ref/S3035822_Covered_sm_lift_hg38.bed";
		
		//String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/vcrome2.1.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/halo_ver2_Regions.bed";
		//String mappability = "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";
		//String mappability = "/GLUSTER_DIST/data/Genomes/karkinos/genome/wgEncodeCrgMapabilityAlign36mer.bw";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinostest2/";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/project/gotoh/S1-S7";
		
		//String outdir = "/GLUSTER_DIST/data/users/ueda/test2";
		//String outdir = "/home/spark/todaitoptest/output/tp0";
		//String outdir = "/home/spark/todaitoptest/karkinos";
		
		String outdir = "/home/spark/todaitoptest/To_Xcoo/test5";
		//String outdir = "/home/spark/todaitoptest/testdata0308/realign/tp0_re";
		String propfile = "/home/spark/todaitoptest/ref/karkinos.properties";
		
		String definedList = "/home/spark/todaitoptest/ref/hotspot.txt";
		
		//
		List<String> l = new ArrayList<String>();
		//add(l,"-n",normalbamf);
		//add(l,"-md",middelfile);
		add(l,"-n",normalbamf);
		add(l,"-t",tumorbamf);
		//add(l,"-chr","chr3");
		//add(l,"-startend","2");
		add(l,"-r",twobitref);
		add(l,"-snp",dbSNP);
		add(l,"-ct",targetRegion);
		//add(l,"-mp",mappability);
		//add(l,"-id","CRC_107_T-CRC_107_N");
		add(l,"-id","TOPTest");
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
