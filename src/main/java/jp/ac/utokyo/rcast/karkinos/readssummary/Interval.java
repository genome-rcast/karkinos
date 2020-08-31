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
package jp.ac.utokyo.rcast.karkinos.readssummary;

public class Interval implements java.io.Serializable{
	
	
	public String getGeneSymbol() {
		return geneSymbol;
	}

	public Interval(String chr,int pos,int depth){
		this.chr = chr;
		this.start = pos;
		this.end = pos;
		this.depth = depth;
	}
	
	public Interval(String chr,int start,int end,String refseqid,String geneSymbol){
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.refseqid = refseqid;
		this.geneSymbol = geneSymbol;
	}
	
	String refseqid = null;
	public int getExonidx() {
		return exonidx;
	}

	String geneSymbol= null;
	int exonidx = 0;
	
	public String getRefseqid() {
		return refseqid;
	}

	public boolean extendInterval(final String chr, final int end, final int depth) {
		if (!this.chr.equals(chr)) {
			return false;
		}
		if (end <= this.end) {
			return false;
		}
		if (this.depth != depth) {
			return false;
		}
		this.end = end;
		return true;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public int getDepth() {
		return depth;
	}

	String chr;
	int start;
	int end;
	int depth;
	public boolean contain(String chr2, int pos) {
		
		boolean bool = chr.equals(chr2)&& (pos>=start && pos<=end);
		return bool;
	}
	
	
}
