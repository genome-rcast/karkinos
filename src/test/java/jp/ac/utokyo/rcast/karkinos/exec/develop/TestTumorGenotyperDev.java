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
package jp.ac.utokyo.rcast.karkinos.exec.develop;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.Option;

import jp.ac.utokyo.rcast.karkinos.exec.TumorGenotyperReanalysis;

public class TestTumorGenotyperDev {
	
	public static void main(String[] arg){
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/sampleN_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/sampleT_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/SNVTest2/C42_tumor_genome.bam";
		
//		String normalbamf = "/data/solexa_data/run_result_hiseq/120328_SN1095_0106_AD0RY3ACXX/Aligned/Project_D0RY3ACXX" +
//				"/Sample_C65_blood_20120326/" +
//				"120328_SN1095_0106_AD0RY3ACXX_C65_blood_20120326_genome.bam";
		
		//String middelfile="/GLUSTER_DIST/data/users/ueda/testframework/testin/sobj";
//		String middelfile="/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_ncc/HCC31T-HCC31N/HCC31T-HCC31N/sobj";
//		String middelfile="/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/THCC_61T/sobj";

//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/C39_tumor_genome.bam";
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/C39_normal_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/project/Koso/20120704_ver0.1/C42recalib_test/C42_normal_genome_recal.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/project/Koso/20120704_ver0.1/C42recalib_test/C42_tumor_genome_recal.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_61T-THCC_61N/normal/THCC_61T-THCC_61N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_61T-THCC_61N/tumor/THCC_61T-THCC_61N_tumor_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_109T-THCC_109N/normal/THCC_109T-THCC_109N_normal_genome.bam";
		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_124T-THCC_124N/tumor/THCC_124T-THCC_124N_tumor_genome.bam";
		String middelfile="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_rcast/THCC_124T-THCC_124N/THCC_124T-THCC_124N/sobj";
		
		
		//String normalbamf = "/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_ncc/HCC31T-HCC31N/normal/HCC31T-HCC31N_normal_genome.bam";
		//String tumorbamf = "/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_ncc/HCC31T-HCC31N/tumor/HCC31T-HCC31N_tumor_genome.bam";
		
		
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		//String dbSNP = "/GLUSTER_DIST/data/users/ueda/SNVtest/hg19_ALL.sites.2012_02.txt";
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_snp132.txt";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_ALL.sites.2012_02.txt";
		
		String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		String mappability = "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";
		
		//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinostest2/";
		String outdir = "/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/THCC_124T-THCC_124N";
		
		//String outdir = "/GLUSTER_DIST/data/users/ueda/KarkinosTest3";
		String propfile = "/GLUSTER_DIST/data/users/ueda/SNVtest/karkinos.properties";
		
		//
		List<String> l = new ArrayList<String>();
		//add(l,"-n",normalbamf);
		add(l,"-md",middelfile);
		add(l,"-t",tumorbamf);
		//add(l,"-chr","1");
		add(l,"-r",twobitref);
		add(l,"-snp",dbSNP);
		add(l,"-ct",targetRegion);
		add(l,"-mp",mappability);
		add(l,"-id","THCC_124T-THCC_124N");
		add(l,"-o",outdir);
		add(l,"-prop",propfile);
		add(l,"-g1000",g1000);
		add(l,"-g1000freq","0.01");
		add(l,"-rg","/GLUSTER_DIST/data/Genomes/ucscgene_hg19/totoki_hg19/refFlat_120402.txt");
		add(l,"-cosmic","/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf");
		add(l,"-exonSNP","/usr/local/karkinos/karkinos/genome/exomeSNP.vcf");
		
		//add(l,"-rs","/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/total_reads_stats_n.txt," +
		//		"/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/total_reads_stats_t.txt");
		
		String[] ar = l.toArray(new String[l.size()]);
		TumorGenotyperDev.main(ar);
		//TumorGenotyper.main(ar);
		
	}

	private static void add(List<String> l, String s1, String s2) {
		l.add(s1);
		l.add(s2);
	}

}
