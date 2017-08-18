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
package jp.ac.utokyo.rcast.karkinos.summary;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;

public class DataWriter {

	public static void writeTr(DataReader dr, XSSFSheet sheet,
			List<String> sampleL) {

		List<String> title = new ArrayList<String>();
		title.add("sample");
		title.add("tumor rate");
		title.add("s.d.");
		title.add("source");
		title.add("number of hetro snp used");
		title.add("correl");
		title.add("ploidy");
		// title.addAll(sampleL);
		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);

		int col = 1;
		for (String sid : sampleL) {
			Row row1 = sheet.createRow(col);
			List data1 = new ArrayList();
			Filebean fb = dr.getBean(sid);

			TCBean tcbean = fb.getTCBean();
			float tr = fb.getTr();
			float sd = fb.getTCBean().sd;
			int nosnp = fb.getTCBean().nosnp;
			float correl = fb.getTCBean().correl;
			float ploidy = fb.getTCBean().ploidy;
			if (tr > 1) {
				tr = 0;
			}
			data1.add(sid);
			data1.add(tr);
			data1.add(sd);
			data1.add("n=" + fb.getTCBean().takefrom);
			data1.add(nosnp);
			data1.add(correl);
			data1.add(ploidy);
			setRow(data1, row1);
			col++;
		}

	}

	public static void writeCNV(DataReader dr, XSSFSheet sheet,
			List<String> sampleL, int flg) {

		List<String> title = new ArrayList<String>();
		title.add("sample id");
		title.add("copy number");
		title.add("chrom");
		title.add("start");
		title.add("end");

		if (flg == 1) {
			title.add("number of hetro SNP");
			title.add("AAF");
			title.add("BAF");
		}

		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);

		int cnt = 1;
		for (String sid : sampleL) {
			Filebean fb = dr.getBean(sid);
			for (CNAInterval cintvl : fb.cnalist) {

				if (flg == 1) {
					if (!cintvl.type.contains("Depth")) {
						continue;
					}
				}
				if (flg == 2) {
					if (!cintvl.type.contains("allelic")) {
						continue;
					}
				}
				if (flg == 3) {
					boolean b = cintvl.type.contains("HD")
							|| cintvl.type.contains("Amp");
					if (!b) {
						continue;
					}
				}

				Row row = sheet.createRow(cnt);
				List data1 = new ArrayList();
				data1.add(sid);
				data1.add(cintvl.cn);
				data1.add(cintvl.chr);
				data1.add(cintvl.start);
				data1.add(cintvl.end);
				if (flg == 1) {
					data1.add(cintvl.noSNP);
					//data1.add(cintvl.varidated);
					data1.add(cintvl.aaf);
					data1.add(cintvl.baf);
				}
				setRow(data1, row);
				cnt++;
			}

		}

	}

	private static void setRow(List sl, Row row) {

		int n = 0;
		for (Object o : sl) {
			Cell c = row.createCell(n);
			if (o instanceof String) {
				c.setCellValue((String) o);
			} else if (o instanceof Double) {
				double d = (Double) o;
				if (d == Double.NEGATIVE_INFINITY) {
					c.setCellValue("n.v.");
				} else {
					c.setCellValue(d);
				}

			} else if (o instanceof Integer) {

				c.setCellValue((Integer) o);

			} else if (o instanceof Long) {

				c.setCellValue((Long) o);

			} else if (o instanceof Float) {

				c.setCellValue((Float) o);

			} else {

				c.setCellValue(String.valueOf(o));
			}
			n++;
		}

	}

	public static void writeReadsStats(DataReader dr, XSSFSheet sheet,
			List<String> sampleL) {
		List<String> title = new ArrayList<String>();
		title.add("type");
		//String anysample = "";
		
		Set<String> keyset = new LinkedHashSet<String>();
		
		for (String s : sampleL) {
			title.add(s + ":normal");
			title.add(s + ":tumor");
			Filebean fb = dr.getBean(s);
			keyset.addAll(fb.datamap.keySet());			
		}
		
		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);

		//Filebean fb = dr.getBean(anysample);
		//Iterator<String> ite = fb.datamap.keySet().iterator();
		
		int col = 1;
		for(String key:keyset) {
			//
			//String key = ite.next();
			List data = new ArrayList();
			data.add(key);
			for (String s : sampleL) {
				Filebean fba = dr.getBean(s);
				String[] adata = fba.datamap.get(key);
				if (adata != null) {
					data.add(toNumber(adata[1]));
					data.add(toNumber(adata[2]));
				}else{
					data.add(0);
					data.add(0);
				}
			}
			Row row = sheet.createRow(col);
			setRow(data, row);
			col++;
		}

	}

	private static Object toNumber(String s) {

		s = s.replaceAll("%", "");
		s = s.replaceAll(",", "");
		//
		try {
			return Long.parseLong(s);
		} catch (Exception ex) {
		}
		try {
			return Float.parseFloat(s);
		} catch (Exception ex) {
		}

		return s;
	}

	public static void writeCNVCBAL(DataReader dr, XSSFSheet sheet,
			List<String> sampleL, ChromBand cband, CellStyle[] csa, boolean high)
			throws SQLException {

		List<String> title = new ArrayList<String>();
		title.add("chr");
		title.add("start");
		title.add("end");
		title.add("name");
		title.add("gieStain");
		title.add("key");
		title.add("mean");
		title.add("sd");			
		title.addAll(sampleL);

		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);
		int col = 1;
		for (CNAInterval band : cband.getList()) {

			List data = new ArrayList();
			data.add(band.chr);
			data.add(band.start);
			data.add(band.end);
			data.add(band.name);
			data.add(band.getGieStain());
			data.add(band.chr + ":" + band.name);
			SummaryStatistics ss = new SummaryStatistics();
		
			for (String sid : sampleL) {

				//
				Filebean fb = dr.getBean(sid);
				float intersectval = 1;
				if(high){
					intersectval = intersect(band, fb ,high);
				}else{
					intersectval = intersect(band, fb ,high);
				}
				
				if(intersectval<=1){
					ss.addValue(intersectval);
				}
				if(intersectval>=1){
					ss.addValue(intersectval);
				}

			}
			data.add(ss.getMean());
			data.add(ss.getStandardDeviation());
			
			for (String sid : sampleL) {

				//
				Filebean fb = dr.getBean(sid);
				float intersectval = 1;
				if(high){
					intersectval = intersect(band, fb ,high);
				}else{
					intersectval = intersect(band, fb ,high);
				}
				
				data.add(intersectval);

			}
			
			Row row = sheet.createRow(col);
			setRow(data, row);
			setCol(data, row, csa);
			col++;
		}

	}
	
	
	private static float intersect(CNAInterval band, Filebean fb,boolean high) {
		
		List<CNAInterval> iclist = new ArrayList<CNAInterval>();

		for (CNAInterval cna : fb.cnalist) {

			if (intersect(cna, band)) {
				
				
				
				boolean b = cna.type.contains("HD")
				|| cna.type.contains("Amp");
				
				if(!b){
					iclist.add(cna);
				}
				
				
				// return cna.cn;
			}

		}
		if (iclist.size() == 0) {
			return 1;
		} else {
			return duality(iclist, band,high);
		}
	}

	private static float duality(List<CNAInterval> iclist, CNAInterval band, boolean high) {
		double area = 0;
		int bandlen = band.end - band.start;
		double areaorg = bandlen;
		for (CNAInterval ci : iclist) {

			//
			int len = ci.end - ci.start;
			if (len >= bandlen)
				len = bandlen;
			float freq =1;
			//
			if(high){
				freq = ci.aaf;
			}else{
				freq = ci.baf;
			}			
			if (freq < 1) {
				area = area - ((1 - freq) * len);
			} else if(freq > 1) {
				area = area + ((freq- 1) * len);
			}

		}
		double a = areaorg + area;
		float cn = (float) (a / bandlen);		
		return cn;
	}
	
	private static float duality(List<CNAInterval> iclist, CNAInterval band) {
		double area = 0;
		int bandlen = band.end - band.start;
		double areaorg = bandlen;
		for (CNAInterval ci : iclist) {

			//
			int len = ci.end - ci.start;
			if (len >= bandlen)
				len = bandlen;
			float cif = ci.cn;
			if (cif < 1) {
				area = area - ((1 - cif) * len);
			} else if(cif > 1) {
				area = area + ((cif - 1) * len);
			}

		}
		double a = areaorg + area;
		float cn = (float) (a / bandlen);		
		return cn;
	}

	public static void writeCNVCB(DataReader dr, XSSFSheet sheet,
			List<String> sampleL, ChromBand cband, CellStyle[] csa, boolean focal)
			throws SQLException {

		List<String> title = new ArrayList<String>();
		title.add("chr");
		title.add("start");
		title.add("end");
		title.add("name");
		title.add("gieStain");
		title.add("key");
		title.add("highmean");
		title.add("sd");
		title.add("lowmean");
		title.add("sd");		
		title.addAll(sampleL);

		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);
		int col = 1;
		for (CNAInterval band : cband.getList()) {

			List data = new ArrayList();
			data.add(band.chr);
			data.add(band.start);
			data.add(band.end);
			data.add(band.name);
			data.add(band.getGieStain());
			data.add(band.chr + ":" + band.name);
			SummaryStatistics sshigh = new SummaryStatistics();
			SummaryStatistics sslow = new SummaryStatistics();
			for (String sid : sampleL) {

				//
				Filebean fb = dr.getBean(sid);
				float intersectval = intersect(band, fb,focal,fb.getTCBean().ploidy);
				if(intersectval<=2){
					sslow.addValue(intersectval);
				}
				if(intersectval>=2){
					sshigh.addValue(intersectval);
				}

			}
			data.add(sshigh.getMean());
			data.add(sshigh.getStandardDeviation());
			data.add(sslow.getMean());
			data.add(sslow.getStandardDeviation());
			
			for (String sid : sampleL) {

				//
				Filebean fb = dr.getBean(sid);
				float intersectval = intersect(band, fb,focal,fb.getTCBean().ploidy);
				data.add(intersectval);

			}
			
			Row row = sheet.createRow(col);
			setRow(data, row);
			setCol(data, row, csa);
			col++;
		}

	}

	private static void setCol(List data, Row row, CellStyle[] csa) {

		int rowno = 0;
		for (Object obj : data) {

			//
			if (obj instanceof Float) {

				float f = (Float) obj;
				if (f > 2) {
					row.getCell(rowno).setCellStyle(csa[0]);
				}
				if (f < 2) {
					row.getCell(rowno).setCellStyle(csa[1]);
				}

			}
			rowno++;
		}
	}

	private static float intersect(CNAInterval band, Filebean fb, boolean focal,float ploidy) {

		List<CNAInterval> iclist = new ArrayList<CNAInterval>();

		for (CNAInterval cna : fb.cnalist) {

			if (intersect(cna, band)) {
				
				
				boolean b = cna.type.contains("HD")
				|| cna.type.contains("Amp");
				
				if(focal){
					
	
					if(b){
						iclist.add(cna);
					}
					
				}else{
					if(!b){
						iclist.add(cna);
					}
				}
				
				// return cna.cn;
			}

		}
		if (iclist.size() == 0) {
			return 2;
		} else {
			return duality(iclist, band,focal,ploidy);
		}
	}

	private static float duality(List<CNAInterval> iclist, CNAInterval band,boolean focal, float ploidy) {

		double area = 0;
		int bandlen = band.end - band.start;
		double areaorg = bandlen * ploidy;
		for (CNAInterval ci : iclist) {

			//
			int len = ci.end - ci.start;
			if (len >= bandlen)
				len = bandlen;
			float cif = ci.cn;
			if (cif < ploidy) {
				area = area - ((ploidy - cif) * len);
			} else {
				area = area + ((cif - ploidy) * len);
			}

		}
		double a = areaorg + area;
		float cn = (float) (a / bandlen);
		if(ploidy==4){
			cn = cn/2;
		}
		if(cn<0)cn=0;
		return cn;
	}

	private static boolean intersect(CNAInterval cna, CNAInterval band) {

		if (!cna.chr.equals(band.chr))
			return false;

		return cna.start < band.end && cna.end > band.start;
	}

}
