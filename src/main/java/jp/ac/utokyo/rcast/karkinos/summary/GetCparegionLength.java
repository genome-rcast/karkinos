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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.TreeMap;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;

public class GetCparegionLength {

	public static long loadTargetFromCapregionFile(File f) throws IOException {

		BufferedReader br = new BufferedReader(new InputStreamReader(
				new FileInputStream(f)));
		int totalcnt = 0;
		long totallen = 0;
		try {

			for (;;) {
				String line = br.readLine();
				if (line == null)
					break;
				totalcnt++;
				String[] sa = line.split("\t");
				String chr = sa[0];
				int start = Integer.parseInt(sa[1]);
				int end = Integer.parseInt(sa[2]);
				int length = end - start;
				totallen = totallen + Math.abs(length);
			}

		} finally {
			br.close();
		}
		System.out.println("target rigion with " + totalcnt + " region and "
				+ +totallen + " bp has loaded");
		return totallen;
	}

}
