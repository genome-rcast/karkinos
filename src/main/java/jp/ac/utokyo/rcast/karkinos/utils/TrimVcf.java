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

public class TrimVcf {

	public static void main(String[] arg) {

		TwoBitGenomeReader tgr = new TwoBitGenomeReader(
				new File(
						"/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/hg38.2bit"));
		String targetRegion = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/S3035822_Covered_sm_lift_hg38_padded_sm.bed";

		DataSet dataset;

		String dir = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/1000genome/GRCh38_positions";
		
		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/1000genome/1000g_phase3_ontarget.vcf";

		
		
		
		try {

			
			File filew = new File(out);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(filew)));

			
			dataset = new DataSet(tgr.readIndex());
			dataset.loadTargetBed(targetRegion, tgr);

			List<String> l = new ArrayList<String>();
			for (int n = 1; n <= 24; n++) {
				for (File f : new File(dir).listFiles()) {

					if (f.getName().contains("tbi"))
						continue;
					if (f.getName().contains("wget"))
						continue;
					
					String ss = n+"";
					if(n==23)ss="X";
					if(n==24)ss="Y";
					
					if (f.getName().contains("chr"+ss+".")){
						l.add(f.getName());
						System.out.println(f.getName());
						break;
				    }

				}
			}
			
			int cnt =0;
			int writecnt = 0;
			
		   int fct =0;
			
			for (String s : l) {
				fct++;
				

				
				System.out.println(s);
				//
				InputStream is = null;
				Reader r = null;
				BufferedReader br = null;
				try {
					is = new GZIPInputStream(new BufferedInputStream(new FileInputStream(dir+"/"+s)));
					r  = new InputStreamReader(is, "MS932");
					br = new BufferedReader(r);
					for (;;) {
						String text = br.readLine();	//改行コードは含まれない
						if (text == null) break;
						
						
						//System.out.println(text);
						if(text.startsWith("#")){
							
							if(fct==1){
								pw.write(text+"\n");
							}
							continue;
							
						}
						String sa[] = text.split("\t");
						String chr = "chr"+sa[0];
						int pos = Integer.parseInt(sa[1]);
						//
						cnt++;
						if(dataset.getCh().getCapInterval(chr, pos)!=null){
							
							pw.write(text+"\n");
							 writecnt++;
													
						}
						

						if(cnt%1000==0){
							System.out.println(fct+"  "+writecnt +"/" +cnt +" " + chr +" "+pos);
							
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



			}

			pw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
