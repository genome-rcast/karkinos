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

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import jp.ac.utokyo.rcast.karkinos.utils.GenotypeKeyUtils;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.apache.commons.math.MathException;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFDataFormat;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

public class SummaryDB {

	Connection memorycon = null;

	PreparedStatement mps = null;

	public SummaryDB(File csv) {

		loadCsv(csv);

	}

	public void loadRef(File ref, String tablename) {
		System.out.println("start loadng");

		try {
			PreparedStatement crs = memorycon.prepareStatement("drop table "
					+ tablename);
			crs.execute();
		} catch (Exception ex) {
		}

		try {

			String sqlcreate = "create table " + tablename
					+ " as select * from csvread('" + ref.getAbsolutePath()
					+ "')";
			PreparedStatement crs = memorycon.prepareStatement(sqlcreate);
			crs.execute();

			System.out.println("end loadng csv");

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void loadCsv(File fin) {

		System.out.println("start loadng");

		try {
			PreparedStatement crs = memorycon
					.prepareStatement("drop table csvdata");
			crs.execute();
		} catch (Exception ex) {
		}

		try {

			Class.forName("org.h2.Driver");
			memorycon = DriverManager.getConnection("jdbc:h2:mem:mydb", "sa",
					"");
			String sqlcreate = "create table csvdata as select * from csvread('"
					+ fin.getAbsolutePath() + "')";
			PreparedStatement crs = memorycon.prepareStatement(sqlcreate);
			crs.execute();

			PreparedStatement crs2 = memorycon
					.prepareStatement("create index idx1 on csvdata(geneSymbol)");
			crs2.execute();

			PreparedStatement crs3 = memorycon
					.prepareStatement("create index idx2 on csvdata(sample_ID)");
			crs3.execute();

			PreparedStatement crs4 = memorycon
					.prepareStatement("create index idx3 on csvdata(ExonicFunc)");
			crs4.execute();

			PreparedStatement crs5 = memorycon
					.prepareStatement("create index idx4 on csvdata(Func)");
			crs5.execute();

			System.out.println("end loadng csv");

		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	List<String> sampleL = null;
	Map<String, Double> backGroundMutationRate = new HashMap<String, Double>();

	public void geneStat(DataReader dr, XSSFSheet sheet, boolean bfilter,
			CellStyle[] cs, boolean aaChange, XSSFWorkbook wb, long caplength,
			boolean HDAmp, int minrecurrent,boolean varidateOnly) throws SQLException {

		sampleL = getSAmpleL();
		setBackGroundMutationRate(sampleL, backGroundMutationRate, bfilter,
				aaChange, caplength);

		List<String> title = new ArrayList<String>();

		title.add("chr");
		title.add("gene");
		title.add("cosmic");

		title.add("distinct pos");

		if (aaChange == false) {

			title.add("synonymous");
			title.add("splice site");
			title.add("AA change");

		}
		//
		title.add("pval");
		// title.add("pval2");
		title.add("log_likehood");
		if (aaChange) {
			title.add("gene cds length");
		} else {
			title.add("gene length");
		}
		title.add("num/kb");

		title.add("observed sample");
		title.add("%sample");
		title.add("Indel");
		title.add("total");
		// mutation
		title.add("LJB_PhyloP_Pred");
		title.add("LJB_SIFT_Pred");
		title.add("LJB_PolyPhen2_Pred");
		title.add("LJB_LRT_Pred");
		title.add("LJB_MutationTaster_Pred");
		int dsratidx = title.size();
		int dendidx = title.size() + sampleL.size();
		//
		title.addAll(sampleL);
		title.add("description");

		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);

		String sql = "select chrom,count(*) as n,group_concat(sample_ID) as m ,geneSymbol from csvdata ";

		if (aaChange) {
			sql = sql
					+ " where ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
					+ "and Func in ('exonic','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + getAndString(sql);
			sql = sql + "  Filter2='PASS'  ";
		}
		
		if(varidateOnly){
			
			sql = sql + getAndString(sql);
			sql = sql + " Validated=TRUE  ";			
		}
		
		sql = sql + " group by geneSymbol ,chrom order by n desc , geneSymbol ";
		PreparedStatement crs = memorycon.prepareStatement(sql);

		ResultSet rs = crs.executeQuery();
		List<List> pList = new ArrayList<List>();
		while (rs.next()) {

			
			List data = new ArrayList();
			String samples = rs.getString("m");
			String chrom = rs.getString("chrom");
			String gene = rs.getString("geneSymbol");
			
			int ncnt = rs.getInt("n");
			//System.out.println(gene + "\t"+ ncnt);
			if (aaChange) {
				if (ncnt < minrecurrent
						|| lessthanP(0.02, ncnt, sampleL.size())) {
					continue;
				}
			} else {
				if (ncnt < minrecurrent
						|| lessthanP(0.04, ncnt, sampleL.size())) {
					continue;
				}
			}
			
			// if(!gene.contains("TTN")){
			// continue;
			// }

			String gene2 = gene;
			if (gene2.contains(";")) {
				gene2 = gene2.substring(0, gene2.indexOf(";"));
			}

			data.add(chrom);
			data.add(gene);
			// cosmic
			String cosmicgene = getStr("select cosmic from generef where geneName ='"
					+ gene2 + "'");
			// data.add(notEmpty(cosmicgene));
			if (cosmicgene == null) {
				cosmicgene = "0";
			}
			data.add(!cosmicgene.equals("0"));
			//
			int total = 0;
			int samplecnt = 0;

			// title.add("Indel");
			// title.add("distinct pos");
			//
			// if(aaChange==false){
			//
			// title.add("synonymous");
			// title.add("splice");
			// title.add("AA change");
			//
			// }

			data.add(getCountSQL(aaChange, bfilter,
					"select count(distinct pos) from csvdata where geneSymbol ='"
							+ gene + "'"));

			if (aaChange == false) {
				data.add(getCountSQL(aaChange, bfilter,
						"select count(*) from csvdata where geneSymbol ='"
								+ gene + "'"
								+ " and ExonicFunc ='synonymous SNV' "));

				data.add(getCountSQL(aaChange, bfilter,
						"select count(*) from csvdata where geneSymbol ='"
								+ gene + "'"
								+ " and Func in ('splicing','exonic;splicing')"));

				data.add(getCountSQL(
						aaChange,
						bfilter,
						"select count(*) from csvdata where geneSymbol ='"
								+ gene
								+ "'"
								+ " and Func in ('exonic','exonic;splicing') and "
								+ " ExonicFunc !='synonymous SNV'"));
			}

			// System.out.println("time1=" +
			// (Calendar.getInstance().getTime().getTime() - start.getTime()));
			// gene length,num/kb
			int genelength = getGenelength(gene2, aaChange);

			double pvaltotal = 0;
			List datatmp = new ArrayList();
			int allscnt = 0;

			double likehood = 0;
			//
			// if(gene2.contains("TGIF2LX")){
			// System.out.println("here");
			// }
			for (String sampleid : sampleL) {
				String tmp = "";
				int countMutation = count(samples, sampleid);
				allscnt++;
				if (countMutation > 0) {
					samplecnt++;
					total = total + countMutation;
					String dispStr = getDispStr(sampleid, gene, aaChange,
							bfilter);
					tmp = dispStr;
					datatmp.add(dispStr);

				} else {

					if (HDAmp) {
						String hdamp = checkHDAmp(dr, sampleid, gene2);
						if (hdamp != null && hdamp.length() > 0) {
							countMutation = 1;
							samplecnt++;
							total = total + countMutation;
							String dispStr = hdamp;
							tmp = dispStr;
							datatmp.add(dispStr);
						} else {
							datatmp.add("0");
						}
					} else {
						datatmp.add("0");
					}

				}

				double pvaleach = getPval(sampleid, backGroundMutationRate,
						gene, aaChange, bfilter, genelength, tmp);
				// System.out.println(sampleid+"\t"+pvaleach+"\t"+likehood);
				likehood = likehood + (-1 * Math.log(pvaleach));

			}

			String kbratio = ratio(
					(double) ((double) total / (double) sampleL.size()),
					((double) genelength / (double) 1000));
			//

			double pval = 1;
			try {
				pval = ProbabilityUtils.getPValbySample(backGroundMutationRate,
						genelength, likehood);

				System.out.println("pval= " + pval + "\t" + gene + "\t"
						+ genelength);

			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// System.out.println("time5=" +
			// (Calendar.getInstance().getTime().getTime() - start.getTime()));
			data.add(pval);
			data.add(likehood);
			data.add(genelength);
			data.add(kbratio);

			//
			data.add(samplecnt);
			data.add(percent(samplecnt, allscnt));

			data.add(getCountSQL(aaChange, bfilter,
					"select count(*) from csvdata where geneSymbol ='" + gene2
							+ "'" + " and Indel ='Yes'"));
			data.add(total);

			// ("LJB_PhyloP_pred");
			data.add(getConutStrSQL(gene2, "LJB_PhyloP_Pred"));
			// ("LJB_SIFT_Pred");
			data.add(getConutStrSQL(gene2, "LJB_SIFT_Pred"));
			// ("LJB_PolyPhen2_Pred");
			data.add(getConutStrSQL(gene2, "LJB_PolyPhen2_Pred"));
			// ("LJB_LTR_Pred");
			data.add(getConutStrSQL(gene2, "LJB_LRT_Pred"));
			// ("LRT_MutationTaster_Pred");
			data.add(getConutStrSQL(gene2, "LJB_MutationTaster_Pred"));
			//
			data.addAll(datatmp);

			// ("description");
			String desc = getStr("select description from generef where geneName ='"
					+ gene2 + "'");
			if (desc == null) {
				desc = "";
			}
			data.add(desc);

			String desccosmic = getStr("select cosmicdesc from generef where geneName ='"
					+ gene2 + "'");
			if (desccosmic == null) {
				desccosmic = "";
			}
			data.add(desccosmic);

			pList.add(data);

			// System.out.println("time6=" +
			// (Calendar.getInstance().getTime().getTime() - start.getTime()));
			//

		}
		rs.close();

		sort(pList, aaChange);

		int colidx = 1;
		for (List data : pList) {

			int maxcolor = data.size();
			Row row = sheet.createRow(colidx);
			setRow(data, row, wb);
			setColor(data, row, cs, dsratidx, dendidx);
			colidx++;
		}

	}

	private boolean lessthanP(double thres, int ncnt, int total) {

		return ((double) ncnt / (double) total) < thres;

	}

	private double getPval(String sampleid,
			Map<String, Double> backGroundMutationRate, String gene,
			boolean aaChange, boolean bfilter, int genelen, String tmp)
			throws SQLException {

		double pval = 0;

		String sql = "select count(*) from csvdata ";
		if (aaChange) {
			sql = sql
					+ " where ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
					+ "and Func in ('exonic','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + getAndString(sql);
			sql = sql + "  Filter2='PASS'  ";
		}
		sql = sql + getAndString(sql);
		sql = sql + " geneSymbol='" + gene + "' and sample_ID='" + sampleid
				+ "'";

		//
		PreparedStatement crs = memorycon.prepareStatement(sql);
		ResultSet rs = crs.executeQuery();
		int cnt = 0;
		if (rs.next()) {
			cnt = rs.getInt(1);
		}
		if (tmp.equals("HD") || tmp.equals("Amp")) {
			cnt++;
		}
		double p = backGroundMutationRate.get(sampleid);
		pval = ProbabilityUtils.getPval(p, cnt, genelen);
		return pval;

	}

	private void setBackGroundMutationRate(List<String> sampleL,
			Map<String, Double> backGroundMutationRate, boolean bfilter,
			boolean aaChange, long caplength) throws SQLException {

		// String sql = ("select * from generef");
		// PreparedStatement crs = memorycon.prepareStatement(sql);
		// ResultSet rs = crs.executeQuery();
		// String chr = "";
		//
		// int totalexon = 0;
		//
		// int start=0;
		// int end=0;
		// while (rs.next()) {
		//
		// String gene = rs.getString("geneName");
		// String gene2 = gene;
		// if (gene2.contains(";")) {
		// gene2 = gene2.substring(0, gene2.indexOf(";"));
		// }
		// int genelength = getGenelength(gene2);
		// totalexon = totalexon + genelength;
		//
		// }
		//
		//
		long total = caplength;
		// if (aaChange) {
		// String s =
		// "select sum(cast(genecdslength as int)) from generef where longestisoform =1";
		// PreparedStatement ps3 = memorycon.prepareStatement(s);
		// ResultSet rs3 = ps3.executeQuery();
		// if (rs3.next()) {
		// total = rs3.getInt(1);
		// }
		// }
		System.out.println("AAchange=" + aaChange + "\t total=" + total);

		for (String sampleid : sampleL) {

			// sample_ID='" + sampleid + "'";
			String sql = "select count(*) from csvdata where sample_ID='"
					+ sampleid + "'";
			// if (aaChange) {
			// sql = sql + getAndString(sql);
			// sql = sql
			// + " ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
			// + "and Func in ('exonic','exonic;splicing') ";
			// }
			if (bfilter) {
				sql = sql + getAndString(sql);
				sql = sql + "  Filter2='PASS' ";
			}

			PreparedStatement crs2 = memorycon.prepareStatement(sql);
			ResultSet rs2 = crs2.executeQuery();
			int mutationcnt = 0;
			if (rs2.next()) {
				mutationcnt = rs2.getInt(1);
			}
			double r = (double) ((double) mutationcnt / (double) total);
			System.out.println("sample=" + sampleid + "\t cnt=" + mutationcnt
					+ "\t exonlength =" + total + "\t r=" + r);
			backGroundMutationRate.put(sampleid, r);
		}

	}

	private String checkHDAmp(DataReader dr, String sid, String gene)
			throws SQLException {

		//
		String sql = ("select * from generef where geneName ='" + gene + "'");

		PreparedStatement crs = memorycon.prepareStatement(sql);
		ResultSet rs = crs.executeQuery();
		String chr = "";
		int start = 0;
		int end = 0;

		if (rs.next()) {

			chr = rs.getString("chrom");
			start = rs.getInt("txstart");
			end = rs.getInt("txend");

		}
		//
		Filebean fb = dr.getBean(sid);

		for (CNAInterval cintvl : fb.cnalist) {

			boolean b = cintvl.type.contains("HD")
					|| cintvl.type.contains("Amp");
			if (!b) {
				continue;
			}
			if (!chr.equals(cintvl.chr)) {
				continue;
			}
			if (cintvl.intercect(start, end)) {
				return cintvl.type;
			}

		}
		return "";
	}

	private void sort(List<List> pList, boolean aaChange) {

		if (aaChange) {
			Collections.sort(pList, new MyCompPval(4));
		} else {
			Collections.sort(pList, new MyCompPval(7));
		}

	}

	private String getDispStr(String sampleid, String gene, boolean aaChange,
			boolean bfilter) throws SQLException {

		String sql = "select * from csvdata ";

		if (aaChange) {
			sql = sql
					+ " where ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
					+ "and Func in ('exonic','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + getAndString(sql);
			sql = sql + "  Filter2='PASS' ";
		}
		sql = sql + getAndString(sql);
		sql = sql + " geneSymbol='" + gene + "' and sample_ID='" + sampleid
				+ "'";

		PreparedStatement crs = memorycon.prepareStatement(sql);
		StringBuffer sb = new StringBuffer();

		ResultSet rs = crs.executeQuery();
		int cnt = 0;
		float cn = 2;
		String dmageStr = "";
		while (rs.next()) {
			cnt++;
			float f = rs.getFloat("CopyNumber");
			if (f != 2) {
				cn = f;
			}
			List<String> slist = new ArrayList<String>();
			slist.add(rs.getString("LJB_PhyloP_Pred"));
			slist.add(rs.getString("LJB_SIFT_Pred"));
			slist.add(rs.getString("LJB_PolyPhen2_Pred"));
			slist.add(rs.getString("LJB_LRT_Pred"));
			slist.add(rs.getString("LJB_MutationTaster_Pred"));
			String dmageStrl = analyse(slist);
			if (dmageStr.equals("") || !dmageStr.equals("D")) {
				if (dmageStrl.equals("D")) {
					dmageStr = dmageStrl;
				} else if (dmageStr.equals("")) {
					dmageStr = dmageStrl;
				}
			}
		}
		if (cnt > 0) {
			sb.append(cnt + ";" + cn + ";" + dmageStr);
		} else if (cn != 2) {
			sb.append(cnt + ";" + cn);
		} else {
			sb.append(0);
		}
		return sb.toString();
	}

	private String analyse(List<String> slist) {

		// C or N or D
		for (String s : slist) {
			if (s == null)
				continue;
			if (s.contains("D") || s.contains("C")) {
				return "D";
			}
			if (s.contains("P")) {
				return "P";
			}
			if (s.contains("B")) {
				return "B";
			}
		}

		return "";

	}

	private boolean notEmpty(String s) {

		return s != null && s.trim().length() > 0;
	}

	private String getStr(String sql) throws SQLException {
		PreparedStatement crs = memorycon.prepareStatement(sql);
		ResultSet rs = crs.executeQuery();
		if (rs.next()) {
			return rs.getString(1);
		}
		return null;
	}

	private String getStrs(String sql) throws SQLException {

		StringBuffer sb = new StringBuffer();
		PreparedStatement crs = memorycon.prepareStatement(sql);
		ResultSet rs = crs.executeQuery();
		while (rs.next()) {
			if (sb.length() > 0) {
				sb.append(",");
			}
			sb.append(rs.getString(1));
		}
		return sb.toString();
	}

	private int getGenelength(String gene, boolean aaChange) {

		int ret = 0;
		try {
			String sql = "select genecdslength from generef where geneName = '"
					+ gene + "'";
			if (!aaChange) {
				sql = "select genelength from generef where geneName = '"
						+ gene + "'";
			}
			PreparedStatement crs = memorycon.prepareStatement(sql);
			ResultSet rs = crs.executeQuery();

			int max = 0;
			while (rs.next()) {
				int n = rs.getInt(1);

				if (max < n) {
					max = n;
				}
			}
			ret = max;
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		System.out.println(gene + "\t" + ret);
		if (ret == 0) {
			ret = 1000;
		}
		return ret;
	}

	// deprocate
	// private int getGeneCDSlength(String gene) {
	//
	// int ret = 0;
	// try {
	// String sql =
	// "select cdsstart,cdsend,exonStarts, exonEnds from generef where geneName = '"
	// + gene + "'";
	// PreparedStatement crs = memorycon.prepareStatement(sql);
	// ResultSet rs = crs.executeQuery();
	//
	// int len = 0;
	// while (rs.next()) {
	//
	// int cdsstart = rs.getInt(1);
	// int cdsend = rs.getInt(2);
	// String starts = rs.getString(3);
	// String ends = rs.getString(4);
	// List<Integer> ls = getIl(starts);
	// List<Integer> le = getIl(ends);
	// int n = 0;
	// for (int ee : le) {
	// int es = ls.get(n);
	//
	// //5'utr
	// if((es<cdsstart) && (es<cdsstart)){
	// n++;
	// continue;
	// }else if(es<cdsstart){
	// es = cdsstart;
	// }
	//
	// //3'utr
	// if((ee>cdsend) && (es>cdsend)){
	// break;
	// }else if(ee>cdsend){
	// ee = cdsend;
	// }
	//
	// len = len + (ee - es);
	// n++;
	// }
	// if (len > ret) {
	// ret = len;
	// }
	// }
	// } catch (Exception ex) {
	// ex.printStackTrace();
	// }
	// if (ret == 0) {
	// ret = 1000;
	// }
	// return ret;
	// }

	private List<Integer> getIl(String s) {
		List<Integer> l = new ArrayList<Integer>();
		for (String ss : s.split(",")) {

			int n = 0;
			try {
				n = Integer.parseInt(ss);
			} catch (Exception ex) {

			}
			l.add(n);
		}
		return l;
	}

	private String getConutStrSQL(String gene, String fld) throws SQLException {

		String sql = "select " + fld
				+ " ,count(*) from csvdata where geneSymbol ='" + gene + "'"
				+ " group by " + fld + " order by " + fld;

		PreparedStatement crs = memorycon.prepareStatement(sql);
		ResultSet rs = crs.executeQuery();
		StringBuffer sb = new StringBuffer();
		while (rs.next()) {
			int cnt = rs.getInt(2);
			String fn = rs.getString(1);
			if (fn == null) {
				continue;
			}
			if (sb.length() > 0) {
				sb.append(",");
			}
			sb.append(fn + cnt);
		}
		return sb.toString();
	}

	private String ratio(double d, double dv) {

		double p = (d / dv);
		return format(p);

	}

	private String percent(double d, double dv) {

		double p = (d / dv) * 100;
		return format(p) + "%";

	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	private int getCountSQL(boolean aaChange, boolean bfilter, String sql) {

		int ret = 0;
		if (aaChange) {
			sql = sql
					+ " and ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
					+ "and Func in ('exonic','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + " and Filter2='PASS'  ";
		}

		try {

			PreparedStatement crs = memorycon.prepareStatement(sql);
			ResultSet rs = crs.executeQuery();
			if (rs.next()) {
				ret = rs.getInt(1);
			}
		} catch (SQLException ex) {

		}
		return ret;
	}

	private void setColor(List data, Row row, CellStyle[] cs, int startindex,
			int endindex) {

		int n = startindex;
		for (; n < endindex; n++) {

			Object obj = data.get(n);
			String s = String.valueOf(obj);
			if (s.contains(";")) {
				String[] sa = s.split(";");
				try {
					float cn = Float.parseFloat(sa[1]);
					if (cn > 2) {
						row.getCell(n).setCellStyle(cs[0]);
					} else if (cn < 2) {
						row.getCell(n).setCellStyle(cs[1]);
					} else {
						row.getCell(n).setCellStyle(cs[2]);
					}
				} catch (Exception ex) {

				}

			} else {

				if (s.equals("HD")) {
					row.getCell(n).setCellStyle(cs[1]);
				} else if (s.equals("Amp")) {
					row.getCell(n).setCellStyle(cs[0]);
				} else if (!s.equals("0")) {
					row.getCell(n).setCellStyle(cs[2]);
				}

			}

		}

	}

	public List<String> getSampleL() {
		return sampleL;
	}

	private List<String> getSAmpleL() throws SQLException {
		if (sampleL != null) {
			return sampleL;
		}
		sampleL = new ArrayList<String>();
		String sqlpre = "select distinct(sample_ID) from csvdata order by sample_ID";
		PreparedStatement crspre = memorycon.prepareStatement(sqlpre);
		ResultSet rspre = crspre.executeQuery();
		while (rspre.next()) {
			sampleL.add(rspre.getString(1));
		}
		rspre.close();
		return sampleL;
	}

	public int mutationStat(XSSFSheet sheet, int colcount, boolean all,
			boolean bfilter) throws SQLException {

		int col = substitution_Stat(sheet, colcount, all, bfilter);
		col = rate_Stat(sheet, col + 5, all, false, bfilter);
		col = rate_Stat(sheet, col + 5, all, true, bfilter);

		return col;
	}

	public int substitution_Stat(XSSFSheet sheet, int colcount, boolean all,
			boolean bfilter) throws SQLException {

		sampleL = getSAmpleL();

		List<String> title = new ArrayList<String>();
		title.add("mutation type\\sample");
		title.addAll(sampleL);
		title.add("total");
		title.add("type");

		Row rowtitle = sheet.createRow(colcount);
		setRow(title, rowtitle);
		colcount++;

		String sql = "select ALT,REF,Indel,count(*) as c,concat(REF,'to',ALT) as m, sample_ID from csvdata ";
		if (!all) {
			sql = sql
//					+ " where ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
//					+ "and Func in ('exonic','splicing','exonic;splicing','splicing') ";
			+ " where ExonicFunc !='synonymous SNV'  "
			+ "and Func in ('exonic','splicing','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + getAndString(sql);
			sql = sql + " Filter2='PASS' ";
		}
		sql = sql + " group by m,sample_ID order by m,sample_ID ";
		PreparedStatement crs = memorycon.prepareStatement(sql);

		ResultSet rs = crs.executeQuery();
		Map<String, Integer> m = new LinkedHashMap<String, Integer>();

		while (rs.next()) {

			int cnt = rs.getInt("c");
			String ref = rs.getString("m");

			String Indel = rs.getString("Indel");
			String sample = rs.getString("sample_ID");
			String key = "";
			if (Indel == null) {
				key = GenotypeKeyUtils.aggrigateKeys(ref);
				// System.out.println(dispkey+"\t"+sample+ "\t"+cnt+"\t"+Indel);
			} else {

				boolean del = rs.getString("REF").length() > rs
						.getString("ALT").length();
				if (del) {
					key = "deletion";
				} else {
					key = "insersion";
				}

			}

			String mkey = sample + ":" + key;
			//
			if (m.containsKey(mkey)) {
				cnt = cnt + m.get(mkey);
				m.put(mkey, cnt);
			} else {
				m.put(mkey, cnt);
			}

		}
		rs.close();

		ArrayList<String> al = new ArrayList<String>();
		al.addAll(Arrays.asList(GenotypeKeyUtils.keys1));
		al.add("insersion");
		al.add("deletion");

		Map<String, Integer> totalbysample = new HashMap<String, Integer>();
		for (String skey : al) {
			ArrayList data = new ArrayList();
			String dispkey = skey;
			try {
				dispkey = GenotypeKeyUtils.toDispKey(skey);
			} catch (ArrayIndexOutOfBoundsException ex) {

			}
			// System.out.print(dispkey+"\t");
			data.add(dispkey);
			int typetotal = 0;
			for (String sample_ID : sampleL) {
				String mkey2 = sample_ID + ":" + skey;
				Integer cnt2 = m.get(mkey2);
				if (cnt2 == null) {
					cnt2 = 0;
				}
				// System.out.print(cnt2+"\t");
				data.add(cnt2);
				typetotal = typetotal + cnt2;
				if (totalbysample.containsKey(sample_ID)) {
					int ll = totalbysample.get(sample_ID);
					ll = ll + cnt2;
					totalbysample.put(sample_ID, ll);
				} else {
					totalbysample.put(sample_ID, cnt2);
				}
			}
			data.add(typetotal);
			if (all) {
				data.add("all pass filter");
			} else {
				data.add("AA alternation");
			}
			Row row = sheet.createRow(colcount);
			setRow(data, row);
			colcount++;
		}
		ArrayList data = new ArrayList();
		data.add("total");
		int nettotal = 0;
		for (String sample_ID : sampleL) {
			int nn = totalbysample.get(sample_ID);
			data.add(nn);
			nettotal = nettotal + nn;
		}
		data.add(nettotal);
		Row row = sheet.createRow(colcount);
		setRow(data, row);

		return colcount;

	}

	private String getAndString(String sql) {
		if (sql.contains("where")) {
			return " and ";
		} else {
			return " where ";
		}
	}

	public int rate_Stat(XSSFSheet sheet, int colcount, boolean all,
			boolean adjusted, boolean bfilter) throws SQLException {

		sampleL = getSAmpleL();

		List<String> title = new ArrayList<String>();
		if (!adjusted) {
			title.add("orginal ratio");
		} else {
			title.add("adjusted ratio");
		}
		title.addAll(sampleL);
		title.add("total");
		title.add("type");

		Row rowtitle = sheet.createRow(colcount);
		setRow(title, rowtitle);
		colcount++;

		String sql = "select count(*) as c, sample_ID, round((AlleleFreq),1) as m  from csvdata ";
		if (!adjusted) {
			sql = "select count(*) as c, sample_ID, round((AlleleFreqOrg),1) as m  from csvdata ";
		}
		if (!all) {
			sql = sql
//					+ " where ExonicFunc !='synonymous SNV' and length(ExonicFunc) > 0 "
//					+ "and Func in ('exonic','splicing','exonic;splicing','splicing') ";
			+ " where ExonicFunc !='synonymous SNV'  "
			+ "and Func in ('exonic','splicing','exonic;splicing','splicing') ";
		}
		if (bfilter) {
			sql = sql + getAndString(sql);
			sql = sql + " Filter2='PASS' ";
		}

		sql = sql + "group by m,sample_ID order by m,sample_ID ";
		PreparedStatement crs = memorycon.prepareStatement(sql);

		ResultSet rs = crs.executeQuery();
		Map<String, Integer> m = new LinkedHashMap<String, Integer>();

		while (rs.next()) {

			int cnt = rs.getInt("c");
			String sampleid = rs.getString("sample_ID");
			float dev = rs.getFloat("m");
			int devi = Math.round(dev * 100);
			if (devi >= 100) {
				devi = 90;
			}
			String key = sampleid + ":" + devi;
			if (m.containsKey(key)) {
				int n = m.get(key);
				cnt = cnt + n;
			}
			m.put(key, cnt);

		}
		rs.close();

		for (int n = 0; n < 100; n = n + 10) {

			ArrayList data = new ArrayList();
			data.add(n + " to " + (n + 10));
			int total = 0;
			for (String sample_ID : sampleL) {

				String key = sample_ID + ":" + n;
				Integer cnt = m.get(key);
				if (cnt == null) {
					cnt = 0;
				}
				total = total + cnt;
				data.add(cnt);
			}
			data.add(total);
			if (all) {
				data.add("all pass filter");
			} else {
				data.add("AA alternation");
			}

			Row row = sheet.createRow(colcount);
			setRow(data, row);
			colcount++;

		}
		return colcount;

	}

	private int count(String samples, String sampleid) {

		int n = 0;
		for (String s : samples.split(",")) {

			if (s.equals(sampleid)) {
				n++;
			}

		}
		return n;
	}

	public void close() throws SQLException {
		if (memorycon != null) {
			try {
				memorycon.createStatement().execute("drop table csvdata");

			} catch (Exception ex) {
			}
			;

			try {
				memorycon.createStatement().execute("drop refflat csvdata");

			} catch (Exception ex) {
			}
			;

			try {
				memorycon.createStatement().execute("drop kgxref csvdata");

			} catch (Exception ex) {
			}
			;

			try {
				memorycon.createStatement().execute("drop cosmic csvdata");

			} catch (Exception ex) {
			}
			;
			if (!memorycon.isClosed()) {
				memorycon.close();
			}
		}
		memorycon = null;
	}

	private static void setRow(List sl, Row row, XSSFWorkbook wb) {

		XSSFCellStyle styleDouble = wb.createCellStyle();
		XSSFDataFormat fd = wb.createDataFormat();
		styleDouble.setDataFormat(fd.getFormat("0.00E+00"));

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
					if (d < 0.01) {
						c.setCellStyle(styleDouble);
					}
				}

			} else if (o instanceof Integer) {

				c.setCellValue((Integer) o);

			} else {

				c.setCellValue(String.valueOf(o));
			}
			n++;
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

			} else {

				c.setCellValue(String.valueOf(o));
			}
			n++;
		}

	}

	public void writeHDAMP(DataReader dr, XSSFSheet sheet, ChromBand cband)
			throws SQLException {

		List<String> title = new ArrayList<String>();
		title.add("sample id");
		title.add("copy number");
		title.add("chrom");
		title.add("start");
		title.add("end");
		title.add("cytoband");
		title.add("length");
		title.add("cosmic genes");
		title.add("#cosmic of genes");
		title.add("genes");
		title.add("#of genes");

		Row rowtitle = sheet.createRow(0);
		setRow(title, rowtitle);

		int cnt = 1;
		for (String sid : sampleL) {
			Filebean fb = dr.getBean(sid);
			for (CNAInterval cintvl : fb.cnalist) {

				boolean b = cintvl.type.contains("HD")
						|| cintvl.type.contains("Amp");
				if (!b) {
					continue;
				}

				Row row = sheet.createRow(cnt);
				List data1 = new ArrayList();
				data1.add(sid);
				data1.add(cintvl.cn);
				data1.add(cintvl.chr);
				data1.add(cintvl.start);
				data1.add(cintvl.end);
				data1.add(cband.getBand(cintvl.chr, cintvl.start, cintvl.end));
				data1.add(Math.abs(cintvl.end - cintvl.start));
				String cgenes = getStrs("select geneName from generef where cosmic is not null and chrom ='"
						+ cintvl.chr
						+ "' and txstart < "
						+ cintvl.end
						+ " and " + cintvl.start + " < txend");

				String genes = getStrs("select geneName from generef where chrom ='"
						+ cintvl.chr
						+ "' and txstart < "
						+ cintvl.end
						+ " and " + cintvl.start + " < txend");
				data1.add(cgenes);
				data1.add(cgenes.split(",").length);
				data1.add(genes);
				data1.add(genes.split(",").length);

				setRow(data1, row);
				cnt++;
			}

		}

	}

	private static String subtype(String type, boolean rev) {

		if (!rev) {
			return type.charAt(0) + "" + type.charAt(1) + "";
		} else {
			return rev(type.charAt(1)) + "" + rev(type.charAt(0)) + "";
		}

	}

	private static char rev(char c) {

		if (c == 'A')
			return 'T';
		if (c == 'T')
			return 'A';
		if (c == 'C')
			return 'G';
		if (c == 'G')
			return 'C';
		return 0;
	}

	public void mutationSig(XSSFSheet sheet, String hg) throws SQLException {

		TwoBitGenomeReader tb = new TwoBitGenomeReader(new File(hg));
		tb.setCheckRepeat(false);
		
		sampleL = getSAmpleL();
		int colcount = 0;
		List<String> title = new ArrayList<String>();
		title.add("mutation type\\sample");
		title.addAll(sampleL);

		Row rowtitle = sheet.createRow(colcount);
		setRow(title, rowtitle);

		String sql = "select sample_ID,chrom,pos,ALT,REF,Indel,concat(REF,'to',ALT) as m from csvdata ";
		sql = sql + " where Filter2='PASS' order by sample_ID,chrom,pos";

		PreparedStatement crs = memorycon.prepareStatement(sql);

		ResultSet rs = crs.executeQuery();
		Map<String, Map<String, Integer>> datam = new LinkedHashMap<String, Map<String, Integer>>();
		Set<String> ts = new TreeSet<String>();
		
		while (rs.next()) {

			String ref = rs.getString("m");
			String sample = rs.getString("sample_ID");
			boolean reverse = ref.startsWith("G") || ref.startsWith("A");
			ref = GenotypeKeyUtils.aggrigateKeys(ref);
			
			String indel = rs.getString("Indel");
			if(indel!=null && indel.equals("Yes"))continue;

			String chr = rs.getString("chrom");
			if(!chr.contains("chr")){
				chr = "chr"+chr;
			}
			int pos = rs.getInt("pos");

			String subseq = "";
			try {
				subseq = tb.getGenomeNuc(chr, pos - 1, true) + ""
						+ tb.getGenomeNuc(chr, pos + 1, true);
			} catch (Exception ex) {
				continue;
			}
			subseq = subtype(subseq, reverse);
			String id = ref + "-" + subseq;
			ts.add(id);
			Map<String, Integer> each = null;
			if (datam.containsKey(sample)) {
				each = datam.get(sample);
			} else {
				each = new HashMap<String, Integer>();
				datam.put(sample, each);
			}
			//
			if (each.containsKey(id)) {
				int nn = each.get(id) + 1;
				each.put(id, nn);
			} else {
				each.put(id, 1);
			}

		}
		rs.close();
		
		for(String muid:ts){
			
			colcount++;
			Row row = sheet.createRow(colcount);
			
			List<String> rowdata = new ArrayList<String>();
			rowdata.add(muid);
			
			for(String sample:sampleL){
				
				Map<String, Integer> each = datam.get(sample);
				Integer num = each.get(muid);
				if(num==null){
					num =0;
				}
				//
				rowdata.add(num+"");				
			}		
			
			setRow(rowdata, row);
			
		}

	}

}
