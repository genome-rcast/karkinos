package jp.ac.utokyo.rcast.karkinos.utils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class TrimVcf2 {

	public static void main(String[] arg) {

		TwoBitGenomeReader tgr = new TwoBitGenomeReader(
				new File(
						"/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/hg38.2bit"));
		String targetRegion = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/S3035822_Covered_sm_lift_hg38_padded_sm.bed";

		DataSet dataset;

//		String dir = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/All_20161122.vcf.gz";
//		
//		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/All_20161122.vcf_ontag.txt";

		String dir = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/Galaxy2-[UCSC_Main_on_Human__snp147Common_(genome)].tabular";
		
		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/snp147Common_hg38.txt";
		
		
		try {

			
			File filew = new File(out);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filew)));

			
			dataset = new DataSet(tgr.readIndex());
			dataset.loadTargetBed(targetRegion, tgr);

			InputStream is = null;
			Reader r = null;
			BufferedReader br = null;
			
			int cnt=0,writecnt=0;
			try {
				is = new BufferedInputStream(new FileInputStream(dir));
				//is = new GZIPInputStream(new BufferedInputStream(new FileInputStream(dir)));
				r  = new InputStreamReader(is, "MS932");
				br = new BufferedReader(r);
				for (;;) {
					String text = br.readLine();	//改行コードは含まれない
					if (text == null) break;
					
					
					//System.out.println(text);
					if(text.startsWith("#")){
						
					
						pw.write(text+"\n");
						
						continue;
						
					}
					String sa[] = text.split("\t");
					String chr = sa[1];
					if(!chr.startsWith("chr")){
						chr = "chr"+chr;
					}
					int pos = Integer.parseInt(sa[2]);
					//
					cnt++;
					if(dataset.getCh().getCapInterval(chr, pos)!=null){
						
						pw.write(text+"\n");
						 writecnt++;
												
					}
					

					if(cnt%1000==0){
						System.out.println(writecnt +"/" +cnt +" " + chr +" "+pos);
						
					}						
					
				}
			} catch (Exception e) {
				//throw new RuntimeException(e);
				e.printStackTrace();
			} finally {
				if (br != null) try { br.close(); } catch (IOException e) {}
				if (r  != null) try { r .close(); } catch (IOException e) {}
				if (is != null) try { is.close(); } catch (IOException e) {}
			}
		

			
				

			pw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
