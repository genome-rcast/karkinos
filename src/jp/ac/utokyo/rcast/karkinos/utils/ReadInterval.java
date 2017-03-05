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

import htsjdk.samtools.SAMRecord;

public class ReadInterval {

	String readname;
	int start;
	int end;
	boolean forward;
	int flg;

	public ReadInterval(SAMRecord sam) {
		add(sam);
	}

	public String getReadname() {
		return readname;
	}

	public void setReadname(String readname) {
		this.readname = readname;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public boolean isForward() {
		return forward;
	}

	public void setForward(boolean forward) {
		this.forward = forward;
	}

	public int getFlg() {
		return flg;
	}

	public void setFlg(int flg) {
		this.flg = flg;
	}

	SAMRecord sam1;
	SAMRecord sam2;
	
	
	public boolean isForward(int pos){
		
		//
		if(sam2==null){
			return forward;
		}
		if((pos-start) < (end-pos) ){
			return true;
		}else{
			return false;
		}
	}
	
	public void add(SAMRecord sam) {

		if(sam1==null){
			sam1 = sam;
		}else{
			sam2 = sam;
		}
		
		int start = sam.getAlignmentStart();
		int end = sam.getAlignmentEnd();
		if (end == 0) {
			end = sam.getAlignmentStart() + sam.getReadLength();
		}

		boolean isfirst = false;
		if (this.start == 0) {
			this.start = start;
			isfirst = true;
		}

		if (this.end == 0) {
			this.end = end;
		}

		if (this.start > start) {
			this.start = start;
		}

		if (this.end < end) {
			this.end = end;
		}

		//
		//
		boolean f1 = false;
		if (sam.getFirstOfPairFlag()) {
			if (sam.getReadNegativeStrandFlag()) {
				f1 = false;
			} else {
				f1 = true;
			}
		} else {
			if (sam.getReadNegativeStrandFlag()) {
				f1 = true;
			} else {
				f1 = false;
			}
		}
//		if (!isfirst) {
//			if (forward == f1) {
//				System.out.println("ok");
//			} else {
//				System.out.println("ng");
//				System.out.println(sam1.format());
//				System.out.println(sam2.format());
//				
//			}
//		}
		forward = f1;

	}

}
