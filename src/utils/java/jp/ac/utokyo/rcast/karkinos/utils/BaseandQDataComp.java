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

import jp.ac.utokyo.rcast.karkinos.filter.BaseandQData;

public class BaseandQDataComp implements Comparator<BaseandQData> {

	public int compare(BaseandQData b0, BaseandQData b1) {
		return b0.getQual()-b1.getQual();
	}
	
}