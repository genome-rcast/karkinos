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

import java.util.Comparator;
import java.util.List;

import org.apache.commons.cli.Option;

public class OptionComparator implements Comparator<Option> {


	List<Option> optionList = null;
	public OptionComparator(List<Option> _optionList) {
		optionList = _optionList;
	}

	public int compare(Option o1, Option o2) {
		int idx1 = optionList.lastIndexOf(o1);
		int idx2 = optionList.lastIndexOf(o2);
		return idx1-idx2;
	}

}
