package jp.ac.utokyo.rcast.karkinos.utils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipInputStream;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

public class Trim1000Data {

	public static void main(String[] arg) {

		String in = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/1000g_phase3_ontarget.vcf";

		//String list = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/nm_ids.txt";

		String out = "/GLUSTER_DIST/data/users/ueda/project/TodaiPanel/Grch38/ref_for_panel/1000g_phase3_ontarget_trim.vcf";

		Set<String> s = new HashSet<String>();

		try {


			File filew = new File(out);
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(
					filew)));

			InputStream is = null;
			Reader r = null;
			BufferedReader br = null;

			int cnt = 0, writecnt = 0;

			Map<String,List<String[]>> holder = new HashMap<String,List<String[]>>();
			
			String title ="";
			int cntn = 0;
			try {
				is = new BufferedInputStream(
						new FileInputStream(in));
				// is = new BufferedInputStream(new FileInputStream(in));
				r = new InputStreamReader(is, "MS932");
				br = new BufferedReader(r);
				int lc = 0;
				for (;;) {
					String text = br.readLine(); // 改行コードは含まれない
					if (text == null)
						break;
					if (text.startsWith("#")) {
						
						pw.write(text+"\n");
						continue;
					}

					String ar[] = text.split("\t");
					cntn++;
					if(cntn%100==0){
						System.out.println(cntn);
					}
					
					for(int n= 0;n<=8;n++){
						
						if(n>0){
							pw.write("\t");
						}
						pw.write(ar[n]);
						
					}					
					pw.write("\n");
					
					lc++;

					
				}
				
				
				
				
			} catch (Exception e) {
				// throw new RuntimeException(e);
				e.printStackTrace();
			} finally {
				if (br != null)
					try {
						br.close();
					} catch (IOException e) {
					}
				if (r != null)
					try {
						r.close();
					} catch (IOException e) {
					}
				if (is != null)
					try {
						is.close();
					} catch (IOException e) {
					}
			}

			pw.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


}
