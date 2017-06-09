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

public class TestTumorGenotyper {
	
	
//	opts.addOption("n", "normalBam", true, "normal bam file");
//	opts.addOption("t", "tumorBam", true, "tumor bam file");
//	opts.addOption("r", "reference", true, "2 bit genome reference file");
//	opts.addOption("snp", "dbSNP", true, "dbSNP list from annover sites, chr,pos,ref,alt,freq,id");
//	opts.addOption("ct", "captureTarget", true, "Capture target regions(bed format)");
//	
//	opts.addOption("o", "outdir", true, "output directory");
//	opts.addOption("id", "uniqueid", true, "unique id for this sample");
//	
//	opts.addOption("prop", "property", true, "path to property file( otherwise use ./karkinos.prop)");
//	opts.addOption("mp", "mappability", true, "optional,mappability from ucsc (bw, big wig format)");
	
	public static void main(String[] arg){
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/sampleN_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/sampleT_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/SNVtest/SNVTest2/C42_tumor_genome.bam";
		
//		String normalbamf = "/data/solexa_data/run_result_hiseq/1203125_SN1095_0106_AD0RY3ACXX/Aligned/Project_D0RY3ACXX" +
//				"/Sample_C65_blood_20120326/" +
//				"1203125_SN1095_0106_AD0RY3ACXX_C65_blood_20120326_genome.bam";
		
		//String middelfile="/GLUSTER_DIST/data/users/ueda/testframework/testin/sobj";
		//String middelfile="/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_ncc/HCC31T-HCC31N/HCC31T-HCC31N/sobj";
//		String middelfile="/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/THCC_61T/sobj";

//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/C39_tumor_genome.bam";
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/karkinostest2/C39/C39_normal_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/project/Koso/20120704_ver0.1/C42recalib_test/C42_normal_genome_recal.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/project/Koso/20120704_ver0.1/C42recalib_test/C42_tumor_genome_recal.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_61T-THCC_61N/normal/THCC_61T-THCC_61N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_61T-THCC_61N/tumor/THCC_61T-THCC_61N_tumor_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_109T-THCC_109N/normal/THCC_109T-THCC_109N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_12-7T-THCC_12-7N/tumor/THCC_12-7T-THCC_12-7N_tumor_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/karkinosSim3/recalbam/ref_recal.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/AC/Ac646T-Ac646N/tumor/Ac646T-Ac646N_tumor_genome.bam";
	
		//String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/bam_baylor/HCC-JP-360-T-HCC-JP-360-N_tumor_genome.bam";
		//String normalbamf = "/GLUSTER_DIST/data/users/ueda/testframework/testkarkinos/THCC_124T/THCC_124T-THCC_124N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/ICGC/HCC_exome_ncc/HCC174T-HCC174N/tumor/HCC174T-HCC174N_tumor_genome.bam";
		
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_RCAST/THCC_125T-THCC_125N/tumor/THCC_125T-THCC_125N_tumor_genome.bam";
//		String middelfile= "/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_ncc/HCC174T-HCC174N/HCC174T-HCC174N/sobj";
		
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/HCC_exome_Baylor/HCC-JP-381-T-HCC-JP-381-N/tumor/HCC-JP-381-T-HCC-JP-381-N_tumor_genome.bam";
//		String middelfile="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor/HCC-JP-186-T-HCC-JP-186-N/HCC-JP-186-T-HCC-JP-186-N/sobj";

		//String tumorbamf = "/GLUSTER_DIST/data/users/ueda/karkinosSim/resultsim1200_rgatkbam/tp90_recal.bam";
//		String middelfile="/GLUSTER_DIST/data/users/ueda/ICGC/ICGCkarkinos4.0/HCC_exome_baylor/HCC-JP-381-T-HCC-JP-381-N/HCC-JP-381-T-HCC-JP-381-N/sobj";
		
		//String middelfile="/GLUSTER_DIST/data/users/ueda/karkinosSim3/karkinos/tp10/tp10/sobj";
		
//		String middelfile="/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/AC/karkinos4.1.6/Ac646T-Ac646N/sobj";

//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/summary_cfDNA1/bam/Se-67-tumor25-Se-67_b-DNA_tumor_genome.bam";
//		String middelfile = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/Se-67-tumor25-Se-67_b-DNA/sobj";
		
		//String middelfile="/GLUSTER_DIST/data/users/yamamoto/work/midorikawa/EHCC_karkinos/karkinos4.1.6/HCC-JP-360-T-HCC-JP-360-N/sobj";
		
//		String tumorbamf = "/GLUSTER_DIST/data/users/gotoh/Hiseq/Halo130408_leukemia_re-analysis/recal/S1_realigned_alnRecal.bam";
		//String middelfile="/GLUSTER_DIST/data/users/ueda/project/gotoh/sobj";
//		String middelfile="/GLUSTER_DIST/data/users/gotoh/Hiseq/Halo130408_leukemia_re-analysis/karkinos/S1-S7/sobj";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/bam_baylor/HCC-JP-418-T-HCC-JP-418-N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/HCC/bam_baylor/HCC-JP-418-T-HCC-JP-418-N_tumor_genome.bam";
		
//		String normalbamf = "/GLUSTER_DIST/data/users/ueda/ICGC/NCCWGS/Normal/HX5WGS_N.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/ueda/ICGC/NCCWGS/Tumor/HX5WGS_T.bam";
	
		//String normalbamf = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/synthetic.challenge.set4.normal_m.bam";
		//String tumorbamf = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/synthetic.challenge.set4.tumor_m.bam";
		
		
//		String normalbamf = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_normal_genome.bam";
//		String tumorbamf = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_tumor_genome.bam";
//		String middelfile="/GLUSTER_DIST/data/users/yamamoto/exome/CRC/karkinos2.0.3/CRC_107_T-CRC_107_N/sobj";
		
		String normalbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3N_TDv3_genome.bam";
		String tumorbamf = "/data/users/yamamoto/TodaiPanel/bam/PLC-TK-3TA_TDv3_genome.bam";
		String middelfile = "/GLUSTER_DIST/data/users/ueda/toptest/sobj/TK-3saveobj.obj";
		
		
//		String normalbamf =  "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT12-1-MT12N/normal/MT12-1-MT12N_normal_genome.bam";
//		String tumorbamf = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT12-1-MT12N/tumor/MT12-1-MT12N_tumor_genome.bam";

		
		//String middelfile = "/data/users/ueda/ICGC/NCCWGS/HX5/sobj";
		
		//String middelfile = "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/karkinos4.1.11/MT/MT12-1-MT12N_k2/sobj";
		
		//String middelfile = "/GLUSTER_DIST/data/users/ueda/project/DreamChallenge/challenge4/sobj";
		
		String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
		//String dbSNP = "/GLUSTER_DIST/data/users/ueda/SNVtest/hg19_ALL.sites.2012_02.txt";
		String dbSNP = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_snp132.txt";
		String g1000 = "/GLUSTER_DIST/data/users/ueda/genome/dbSNP132/hg19_ALL.sites.2012_02.txt";
		
		//String targetRegion = "/data/users/ueda/project/TodaiPanel/S3035822_Covered.bed_test";
		String targetRegion = "/data/users/yamamoto/TodaiPanel/target/S3035822_Covered.bed";
		
		//String targetRegion = "/GLUSTER_DIST/data/Genomes/hg19capv4/captureregionv4all.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/vcrome2.1.bed";
		//String targetRegion ="/GLUSTER_DIST/data/Genomes/karkinos/genome/halo_ver2_Regions.bed";
		//String mappability = "/GLUSTER_DIST/data/users/ueda/SNVtest/wgEncodeCrgMapabilityAlign100mer.bw";
		String mappability = "/GLUSTER_DIST/data/Genomes/karkinos/genome/wgEncodeCrgMapabilityAlign36mer.bw";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/karkinostest2/";
		//String outdir = "/GLUSTER_DIST/data/users/ueda/project/gotoh/S1-S7";
		
		//String outdir = "/GLUSTER_DIST/data/users/ueda/test2";
		String outdir = "/GLUSTER_DIST/data/users/ueda/toptest";
		String propfile = "/usr/local/karkinos/karkinos4.0/karkinos.properties";
		
		String definedList = "/GLUSTER_DIST/data/users/ueda/toptest/deftest.txt";
		
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
		add(l,"-mp",mappability);
		//add(l,"-id","CRC_107_T-CRC_107_N");
		add(l,"-id","TK-3");
		add(l,"-o",outdir);
		add(l,"-prop",propfile);
		add(l,"-g1000",g1000);
		add(l,"-g1000freq","0.01");
		//add(l,"-tc","1.0");
		//add(l,"-rg","/GLUSTER_DIST/data/Genomes/ucscgene_hg19/totoki_hg19/refFlat_120402.txt");
		
		add(l,"-rg","/GLUSTER_DIST/data/users/ueda/toptest/top_refflathg19.txt");
		add(l,"-cosmic","/usr/local/karkinos/karkinos/genome/CosmicCompleteExport_v62_291112.vcf");
		add(l,"-exonSNP","/usr/local/karkinos/karkinos/genome/exomeSNP.vcf");
		
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
