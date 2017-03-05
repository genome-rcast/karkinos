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
package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import java.util.Comparator;

import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;

public class MyComp implements Comparator<CopyNumberInterval> {

	public int compare(CopyNumberInterval arg0, CopyNumberInterval arg1) {

		if (arg0.getChr().equals(arg1.getChr())) {
			return arg0.getStart() - arg1.getStart();
		} else {
			return getIndex(arg0.getChr()) - getIndex(arg1.getChr());
		}

	}

	private int getIndex(String chr) {

		chr = chr.replace("chr", "");
		if (isNumber(chr)) {
			return Integer.parseInt(chr);
		} else if (chr.equals("X")) {
			return 101;
		} else if (chr.equals("Y")) {
			return 102;
		} else if (chr.equals("M")) {
			return 103;
		}
		return chr.hashCode();

	}

	private boolean isNumber(String chr) {

		try {

			int n = Integer.parseInt(chr);
			return true;

		} catch (Exception ex) {
			return false;
		}

	}

}
