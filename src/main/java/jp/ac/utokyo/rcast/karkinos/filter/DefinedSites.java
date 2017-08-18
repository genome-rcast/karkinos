package jp.ac.utokyo.rcast.karkinos.filter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.Set;

public class DefinedSites {

	Set<Integer> positions = new LinkedHashSet<Integer>();

	String sites = null;
	String chr = null;
	public DefinedSites(String sites) {
		
		this.sites = sites;
		
	}

	public boolean contains(String chr,int n) {
		
		if(this.chr==null || this.chr.equals(chr)){
			load(sites,chr); 
		}		
		//
		return positions.contains(n);
	}

	
	
	public void load(String sites, String chrom) {

		this.chr = chr;
		//
		File file = new File(sites);
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(file));

			String str = null;
			while ((str = br.readLine())!= null) {
			
				//
				String[] ary = str.split("\t");
				
				String chr = ary[0];
				if(chr.equals(chrom)){
					
					int pos = Integer.parseInt(ary[1]);
					positions.add(pos);
				}				
				
			}
			br.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public boolean contains(int n) {
		//
		return positions.contains(n);
	}

}
