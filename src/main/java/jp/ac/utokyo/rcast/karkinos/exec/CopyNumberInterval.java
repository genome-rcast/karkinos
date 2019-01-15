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

/**
 * The range is 1-based closed range [start,end]
 */
public class CopyNumberInterval implements java.io.Serializable{
	private final String chr;
	private final float aaf, baf;
	private int start, end;
	private float copynumber;
	private boolean allelic = false;
	private boolean hdeletion = false; // Homozygous deletion
	private int noSNP;

	public CopyNumberInterval(final String chr) {
		this(chr, 0, 0);
	}

	public CopyNumberInterval(final String chr, final float aaf, final float baf) {
		this.chr = chr;
		this.aaf = aaf;
		this.baf = baf;
	}

	public void setNoSNP(final int noSNP) {
		this.noSNP = noSNP;
	}

	public float getAaf() {
		return aaf;
	}

	public float getBaf() {
		return baf;
	}

	public boolean isAllelic() {
		return allelic;
	}

	public void setAllelic(final boolean allelic) {
		this.allelic = allelic;
	}

	public boolean isHdeletion() {
		return hdeletion;
	}

	public void setHdeletion(final boolean hdeletion) {
		this.hdeletion = hdeletion;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public void setStart(final int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(final int end) {
		this.end = end;
	}

	public float getCopynumber() {
		return copynumber;
	}

	public int getNoSNP() {
		return noSNP;
	}

	public void setCopynumber(final float copynumber) {
		this.copynumber = copynumber;
	}

	/**
	 * Return the length of Copy Number Interval.
	 */
	public int length() {
		if (start > end) {
			return 0;
		}
		return end - start + 1;
	}
}
