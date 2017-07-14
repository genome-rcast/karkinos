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

import java.io.EOFException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

public class TwoBitGenomeReader {

	public void setCheckNmask(boolean checkNmask) {
		this.checkNmask = checkNmask;
	}

	public void setCheckRepeat(boolean checkRepeat) {
		this.checkRepeat = checkRepeat;
	}

	File twoBitRef = null;
	boolean checkNmask = true;
	boolean checkRepeat = true;

	Map<String, Integer> index = null;

	public TwoBitGenomeReader(File f) {

		twoBitRef = f;
		try {
			index = readIndex();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public Map<String, Integer> readIndex() throws IOException {

		Map<String, Integer> m = new LinkedHashMap<String, Integer>();
		RandomAccessFile raf = new RandomAccessFile(twoBitRef, "r");
		LittleEndian la = new LittleEndian(raf);
		raf.seek(0L);
		// header
		int sig = la.readInt();
		boolean sigOK = (sig == 0x1a412743);
		int version = la.readInt();
		int seqCount = la.readInt();
		int reserved = la.readInt();

		for (int n = 0; n < seqCount; n++) {
			byte namesize = 0;
			try{
				namesize = raf.readByte();
			}catch(EOFException eof){
				continue;
			}
			byte[] ba = new byte[namesize];
			raf.read(ba);
			String name = new String(ba);
			int size = la.readInt();
			m.put(name, size);
			//System.out.println(name+":"+size);
		}
		raf.close();
		return m;

	}
	
	public Map<String, Integer> getReadSizes() throws IOException {

		Map<String, Integer> m = new LinkedHashMap<String, Integer>();
		Map<String, Integer> m2 = new LinkedHashMap<String, Integer>();
		RandomAccessFile raf = new RandomAccessFile(twoBitRef, "r");
		LittleEndian la = new LittleEndian(raf);
		raf.seek(0L);
		// header
		int sig = la.readInt();
		boolean sigOK = (sig == 0x1a412743);
		int version = la.readInt();
		int seqCount = la.readInt();
		int reserved = la.readInt();

		for (int n = 0; n < seqCount; n++) {
			byte namesize = raf.readByte();
			byte[] ba = new byte[namesize];
			raf.read(ba);
			String name = new String(ba);
			int size = la.readInt();
			m.put(name, size);
		}
		Iterator<String> ite = m.keySet().iterator();
		while(ite.hasNext()){
			String chrom = ite.next();
			raf.seek(m.get(chrom));
			int dnaSize = la.readInt();
			m2.put(chrom, dnaSize);			
		}		
		raf.close();
		return m2;

	}
	
	
	
	public static Map<String, Integer> getReadSizes(File f) throws IOException {

		Map<String, Integer> m = new LinkedHashMap<String, Integer>();
		Map<String, Integer> m2 = new LinkedHashMap<String, Integer>();
		RandomAccessFile raf = new RandomAccessFile(f, "r");
		LittleEndian la = new LittleEndian(raf);
		raf.seek(0L);
		// header
		int sig = la.readInt();
		boolean sigOK = (sig == 0x1a412743);
		int version = la.readInt();
		int seqCount = la.readInt();
		int reserved = la.readInt();

		for (int n = 0; n < seqCount; n++) {
			byte namesize = raf.readByte();
			byte[] ba = new byte[namesize];
			raf.read(ba);
			String name = new String(ba);
			int size = la.readInt();
			m.put(name, size);
		}
		Iterator<String> ite = m.keySet().iterator();
		while(ite.hasNext()){
			String chrom = ite.next();
			raf.seek(m.get(chrom));
			int dnaSize = la.readInt();
			m2.put(chrom, dnaSize);			
		}		
		raf.close();
		return m2;

	}
	
	

	private byte genomeBuf[] = null;
	private String lastloadchr = "";
	private PositionBrocks nBlock = null, maskBlock = null;

	
	public void readRefIfNew(String chrom) throws IOException {
		if (!lastloadchr.equals(chrom)) {
			System.out.println(chrom);
			readRef(chrom);
		}		
	}
	
	public boolean isRefExsist(String chrom){
		
		boolean contain = index.containsKey(chrom);
		if(!contain){
			contain = index.containsKey("chr"+chrom);
		}
		return contain;
	}
	
	public String getChromString(int chridxn) {
		
		int n = 1;
		for(Entry<String, Integer> et:index.entrySet()){
			
			if(chridxn==n){
				return et.getKey();
			}
			n++;
		}		
		return null;
	}
	
	public boolean readRef(String chrom) throws IOException {
		
		//System.out.println(chrom);
		RandomAccessFile raf = new RandomAccessFile(twoBitRef, "r");
		LittleEndian la = new LittleEndian(raf);
		if(!isRefExsist(chrom)){
			//should not happen
			genomeBuf=null;
			lastloadchr = chrom;
			return false;
		}		
		raf.seek(getIdx(index,chrom));
		int dnaSize = la.readInt();
		int nBlockCount = la.readInt();
		// nblock
		if (checkNmask) {
			nBlock = new PositionBrocks();
			nBlock.read(nBlockCount, la);
		} else {
			raf.skipBytes(nBlockCount * 8);
		}
		int maskBlockCount = la.readInt();
		// repeat mask
		if (checkRepeat) {
			maskBlock = new PositionBrocks();
			maskBlock.read(maskBlockCount, la);
		} else {
			raf.skipBytes(maskBlockCount * 8);
		}
		int asize = (dnaSize / 4) + 1;
		int reserved = la.readInt();
		genomeBuf = new byte[asize];
		raf.read(genomeBuf);
		raf.close();
		lastloadchr = chrom;
		return true;
	}

	private long getIdx(Map<String, Integer> index, String chrom) {
		
		boolean contain = index.containsKey(chrom);
		if(contain){
			return index.get(chrom);
		}
		Integer i = index.get("chr"+chrom);		
		if(i==null){
			i=0;
		}
		return i;
		
	}

	public char getGenomeNuc(String chrom, int pos, boolean strand)
			throws IOException {

		if (genomeBuf == null || (!lastloadchr.equals(chrom)&&(!chrom.equals("*")))) {
			System.out.println(chrom);
			readRef(chrom);
		}
		return getGenomeNuc(pos, strand);
	}
	
	public float getCGParcent(String chrom, int start, int end) throws IOException{
		
		int total = end-start;
		float cg = 0;
		for(int n= start;n<=end;n++){
			char c = getGenomeNuc(chrom,n,true);
			if((c=='c')||(c=='C')||(c=='g')||(c=='G')){
				cg++;
			}else if(c=='N'){
				cg = cg+0.5f;
			}
		}
		double r = (double)((double)cg/(double)total)*100;
		return (float)r;
	}

	public char getGenomeNuc(int pos, boolean strand) {

		char nuc = _getGenomeNuc(pos, strand);
		//NMask
		if (checkNmask&&nBlock.inBlock(pos)) nuc= 'N';
		//Repeak
		if (checkRepeat&&maskBlock.inBlock(pos)) {
			return Character.toLowerCase(nuc);
		} else {
			return nuc;
		}

	}

	public char _getGenomeNuc(int pos, boolean strand) {

		if (genomeBuf == null)return 0;
		
		try {
			//
			int arypos = (pos - 1) / 4;
			//
			if(arypos>=genomeBuf.length){
				return 0;
			}
			byte byteofpos = genomeBuf[arypos];
			int narrowPos = (pos - 1) % 4;
			// T-00,C-01,A-10,G-11
			return getNuc(byteofpos, narrowPos, strand);

		} catch (Exception ex) {
			ex.printStackTrace();
		}
		return 0;
	}

	public String getGenomicSeq(String chrom, int start, int end, boolean strand)
			throws IOException {

		StringBuffer sb = new StringBuffer();
		try{
			for (int n = start; n <= end; n++)
				sb.append(getGenomeNuc(chrom, n, strand));
		}catch (ArrayIndexOutOfBoundsException ex){}	
		return strand ? sb.toString() : sb.reverse().toString();

	}

	private char getNuc(byte byteofpos, int narrowPos, boolean strand) {

		int nucint = getNucInt(byteofpos, narrowPos);
		if (strand) {
			switch (nucint) {
			case 0:
				return 'T';
			case 1:
				return 'C';
			case 2:
				return 'A';
			case 3:
				return 'G';
			}
		} else {
			switch (nucint) {
			case 0:
				return 'A';
			case 1:
				return 'G';
			case 2:
				return 'T';
			case 3:
				return 'C';
			}

		}
		return 0;
	}

	private int getNucInt(byte byteofpos, int narrowPos) {

		switch (narrowPos) {
		case 0:
			return (byteofpos & 192) >> 6;
		case 1:
			return (byteofpos & 48) >> 4;
		case 2:
			return (byteofpos & 12) >> 2;
		case 3:
			return (byteofpos & 3);
		}
		return 0;
	}

	public static void main(String[] args) {
		//
		File f = new File("/data/Genomes/hg18_mito/hg18.2bit");
		TwoBitGenomeReader inst = new TwoBitGenomeReader(f);
		try {

			String seq = inst.getGenomicSeq("chr1",100000,100550,true);
				
			System.out.println(seq);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	

	

}
