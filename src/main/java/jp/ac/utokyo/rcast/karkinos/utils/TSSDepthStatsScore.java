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
package jp.ac.utokyo.rcast.karkinos.utils;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.transform.FastFourierTransformer;

public class TSSDepthStatsScore extends ReadWriteBase {

	final static int SIZE = 1024;
	public static void main(String[] arg) throws IOException {

		
		String bamn = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/normal/"
			+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_normal_genome.bam";
		
		String bamt = "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/cfDNA/Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS/tumor/"
			+ "Se68_cfDNA_25ng_WGS-Se68_bDNA_25ng_WGS_tumor_genome.bam";
	

		String bandf = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS.txt";
		String out = "/GLUSTER_DIST/data/users/ueda/project/Asada/RefSeqTSS_Score.txt";

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		
		
		SAMFileReader bamrn = getReader(bamn);
		SAMFileReader bamrt = getReader(bamt);
		

		
		// /
		String line = null;
		
		BufferedReader sband = new BufferedReader(new FileReader(bandf));
		Set<String> kset = new HashSet<String>();
		
		while ((line = sband.readLine()) != null) {

			//if(cntn>10)break;
			if (line.startsWith("#")) {
				continue;
			}
			//debug
//			if (!line.contains("BUB1")) {
//				continue;
//			}

			String[] sa = line.split("\t");

			String ucsc = sa[0];
			
			
			boolean strand = sa[2].trim().equals("+");
			String chr = sa[1];
			int pos = Integer.parseInt(sa[5]);
			String stranbd = sa[3];
			int start = pos -  (SIZE-1);
			int end = pos + SIZE;
			String key = chr + "_" + pos;

			if (kset.contains(key)) {
				continue;
			}
			kset.add(key);
			
			
			double[] forC = new double[SIZE*2];
			double[] revC = new double[SIZE*2];
			
			double[] totalC = new double[SIZE*2];
			

			
			CloseableIterator<SAMRecord> iten = bamrn.query(chr, start, end,
					false);
			List<ReadInterval> rin = getList(iten); 
			
	
			CloseableIterator<SAMRecord> itet = bamrt.query(chr, start, end,
					false);
			
		
			
			List<ReadInterval> rit = getList(itet); 
			
			double dt0 = 0;
			double dt1 = 0;
			
			double dn0 = 0;
			double dn1 = 0;
			
			if (strand) {

				for (int n = -(SIZE-1); n <= SIZE; n++) {
					int pos2 = pos + n;
					int[] retn = getD(pos2, rin);
					int[] rett = getD(pos2, rit);
					double r0 = (double)(rett[0]+1)/(double)(retn[0]+1); 
					double r1 = (double)(rett[1]+1)/(double)(retn[1]+1); 
					
					double rall = (double)(rett[0]+rett[1]+1)/(double)(retn[0]+retn[1]+1); 
//					if(r0>1)r0=1;
//					if(r1>1)r1=1;
					
					forC[n+(SIZE-1)] = r0;
					revC[n+(SIZE-1)] = r1;
					
					totalC[n+(SIZE-1)] = rall;
					
					if(Math.abs(n-10)<=10){
						dn0 = dn0 + retn[0];
						dn1 = dn1 + retn[1];
					}
					
					if(Math.abs(n-10)<=10){
						dt0 = dt0 + rett[0];
						dt1 = dt1 + rett[1];
					}
					
				}

			} else {

				for (int n = SIZE; n >= -(SIZE-1); n--) {
					int pos2 = pos + n;
					int[] retn = getD(pos2, rin);
					int[] rett = getD(pos2, rit);
					
					double r0 = (double)(rett[0]+1)/(double)(retn[0]+1); 
					double r1 = (double)(rett[1]+1)/(double)(retn[1]+1); 	
					
					double rall = (double)(rett[0]+rett[1]+1)/(double)(retn[0]+retn[1]+1); 
//					if(r0>1)r0=1;
//					if(r1>1)r1=1;
//					
					forC[n+(SIZE-1)] = r0;
					revC[n+(SIZE-1)] = r1;
					
					totalC[n+(SIZE-1)] = rall;
					
					if(Math.abs(n-10)<=10){
						dn0 = dn0 + retn[0];
						dn1 = dn1 + retn[1];
					}
					
					if(Math.abs(n-10)<=10){
						dt0 = dt0 + rett[0];
						dt1 = dt1 + rett[1];
					}

				}

			}

			//norm
			
		
//			for (int n = 0; n <= 2048; n++) {
//				
//				forCs[n] = sm(forC,n);
//				revCs[n] = sm(revC,n);
//				
//			}
			
			
//			for (int n = 0; n < SIZE*2; n++) {
//				
//				if((n<512) || (n>(1024+512))){
//				    forC[0] = 2;
//				    revC[0] = 2;
//				}else{
//					 forC[n] = 0;
//					 revC[n] = 0;
//				}
				
//				if(n<=1024){
//			     forC[0] = 1-(double)(n)/1024;
//			     revC[0] = 1-(double)(n)/1024;
//				}else{
//				  forC[0] = -1+(double)(n)/1024;
//				  revC[0] = -1+(double)(n)/1024;
//				}
				
//				if((n/10)%2==0){
//				    forC[n] = 1;
//				    revC[n] = 1;
//				}else{
//				    forC[n] = -1;
//				    revC[n] = -1;
//				}
//				double dd = -Math.sin(2*Math.PI*((double)0.5*n/(double)2048));
//			    forC[n] = dd;
//			    revC[n] = dd;
			    //System.out.println(dd);
//			    forC[n] = 1;
//			    revC[n] = 1;
//			}
			
			FastFourierTransformer fft = new FastFourierTransformer();
			
			Complex[] fw = fft.transform2(forC);
			Complex[] rv = fft.transform2(revC);
			
			////
			////
//			int cnt = 0;
//			for(Complex cp:rv){
//				
//				cnt++;
//				double real = cp.getReal();
//				double im = cp.getImaginary();
//				double mag = Math.sqrt((real*real)+(im*im));
//				
//		
//				if(cnt>SIZE*2) break;
//				System.out.println(mag);
//				
//			}
			
//			for(Complex cp:fw){
//				
//				double real = cp.getReal();
//				double im = cp.getImaginary();
//				double mag = Math.sqrt((real*real)+(im*im));
//				System.out.println(mag);
//				
//			}
			String[] scf = sc(fw);
			String[] scr = sc(rv);
			
			//check concaveup
			double dd = (dt0+dt1)/(dn0+dn1);
			//bw.write((dn0/20)+"\t"+(dn1/20)+"\t"+(dt0/20)+"\t"+(dt1/20)+"\t"+dd+"\t"+scf[0]+"\t"+scf[1]+"\t"+scf[2]+"\t"+scr[0]+"\t"+scr[1]+"\t"+scr[2]+"\t"+line);
			
			//regression
			//-512 to 512
			//find a 
			double [] regression = getRegression(totalC);
			
			StringBuffer sb = new StringBuffer();
			sb.append((dn0/20)+"\t"+(dn1/20)+"\t"+(dt0/20)+"\t"+(dt1/20)+"\t"+dd+"\t");
			sb.append(regression[0]+"\t"+regression[1]+"\t"+regression[2]+"\t"+regression[3]+"\t");
			for(int n=0;n<20;n++){
				sb.append(mag(fw[n])+"\t");
			}
			for(int n=0;n<20;n++){
				sb.append(mag(rv[n])+"\t");
			}
			sb.append(getGS(line)+"\t"+line);
			bw.write(sb.toString()+"\n");
			System.out.println(sb.toString());
			//break;
		}
		
		bw.close();

	}
	
	private static List<ReadInterval> getList(CloseableIterator<SAMRecord> itet) {
		
		
		Map<String,ReadInterval> map = new HashMap<String,ReadInterval>();
			
		while (itet.hasNext()) {

			SAMRecord sam = itet.next();
			if(sam.getReadUnmappedFlag())continue;
			if(!sam.getProperPairFlag())continue;
			//
			ReadInterval ri = null;
			if(map.containsKey(sam.getReadName())){
				
				//
				ri = map.get(sam.getReadName());
				ri.add(sam);
				
			}else{
				
				ri = new ReadInterval(sam);
				map.put(sam.getReadName(),ri);
				
			}
			
		}
		itet.close();
		//
		//
		List<ReadInterval> list = new ArrayList<ReadInterval>();
		Iterator<String> ite = map.keySet().iterator();
		while(ite.hasNext()){
			
			//
			list.add(map.get(ite.next()));
			
		}		
		return list;
		
	}

	private static double[] getRegression(double[] totalC) {
		
				
		//1.regression by 2nd degree polynomial
		int n = 0;
		int m=0;
		double X = 0;
		double X2 = 0;
		double X3 = 0;
		double X4 = 0;
		
		double Y = 0;
		double XY=0;
		double X2Y = 0;
		
		
		for(double d:totalC){
			
			n++;
			
			double x0 = (n-SIZE);
			if(x0>(SIZE/5))continue;
			if(x0<-(SIZE/5))continue;
			
			m++;
			double y0 = d;
			
			//System.out.println(x0+"\t"+y0);
					
			double x2 = Math.pow(x0, 2);
			double x3 = Math.pow(x0, 3);
			double x4 = x2 * x2;

			X += x0;
			X2 += x2;
			X3 += x3;
			X4 += x4;
			

			X2 += x2;
			X4 += x4;
			
			Y += y0;
			XY += x0 * y0;
			X2Y += x2 * y0;
			
		}
		
		// resolve equation
		double[][] vals = { { n, X, X2 }, { X, X2, X3 }, { X2, X3, X4 } };
		double[] rhs = { Y, XY, X2Y };

		RealMatrix rm = new Array2DRowRealMatrix(vals);

		double a1 = 0;
		double b1 = 0;
		double c1 = 0;


		try {
			DecompositionSolver solver = new LUDecompositionImpl(rm).getSolver();
			RealVector b = new ArrayRealVector(rhs);
			RealVector x = solver.solve(b);
			a1 = x.getEntry(2);
			b1 = x.getEntry(1);
			c1 = x.getEntry(0);
	
		} catch (Exception ex) {
			//ex.printStackTrace();
		}
		
		double red = 0;
		n=0;
		for(double d:totalC){
			n++;
			double x0 = (n-SIZE);
			if(x0>(SIZE/5))continue;
			if(x0<-(SIZE/5))continue;
			double y0 = d;			
			double yc = a1*(x0*x0)+b1*(x0)+c1;
			
			red = red + Math.pow((y0-yc),2);
			
		}
		if(red>50){
			a1=0;
			b1=0;
			c1=0;
			red =0;
		}
		return new double[]{a1,b1,c1,red};
	}

	private static String getGS(String line) {
		
		try{
			return line.substring(line.lastIndexOf("(")+1,line.lastIndexOf(")"));
		}catch(Exception ex){
			
		}
		return "";
	}

	private static String[] sc(Complex[] cpa) {

		double total = 0;
		for(Complex c: cpa){
			double d = mag(c);
			total = total + d;
		}		
		
		String[] ret = new String[3];
		double s1 = mag(cpa[0])/total;
		
		
		ret[0]=s1+"";
		ret[1]=(mag(cpa[2])/total+mag(cpa[3])/total+mag(cpa[4])/total+mag(cpa[5])/total+mag(cpa[6])/total)+"";
		ret[2]=(mag(cpa[21])/total+mag(cpa[22])/total+mag(cpa[23])/total+mag(cpa[24])/total+mag(cpa[25])/total)+"";
		return ret;
		
	}



	private static double mag(Complex cp) {
		double real = cp.getReal();
		double im = cp.getImaginary();
		double mag = Math.sqrt((real*real)+(im*im));
		return mag;
	}

	private static int[] getD(int pos2, List<ReadInterval> rin) {

		int f = 0;
		int r = 0;
		for (ReadInterval ri : rin) {


			if (ri.getStart() > pos2)
				continue;
			if (ri.getEnd()  < pos2)
				continue;
			
			if(ri.isForward(pos2)){
				f++;
			}else{
				r++;
			}

		}		
		//
		return new int[] { f, r };

	}
	
	
}
