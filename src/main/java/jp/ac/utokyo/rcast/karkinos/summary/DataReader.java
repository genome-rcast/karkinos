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

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import au.com.bytecode.opencsv.CSVReader;

public class DataReader {

	List<Filebean> list;

	public DataReader(List<Filebean> list) throws IOException {

		this.list = list;
		for (Filebean fb : list) {

			float tr = 0f;
			float sd = 0f;
			int nohsnp = 0;
			float correl = 0f;

			CSVReader brcvs = new CSVReader(new FileReader(fb.tdata), '\t');
			String[] data = null;
			List<CNAInterval> l = new ArrayList<CNAInterval>();

			Map<String, String[]> readstats = new LinkedHashMap<String, String[]>();

			while ((data = brcvs.readNext()) != null) {
				//

				if (data != null
						&& (data[0].equals("n=1") || data[0].equals("n=3") || data[0]
								.equals("somatic"))) {

					String name = data[0];
					tr = toFloat(data[1]);
					sd = toFloat(data[2]);
					nohsnp = toInt(data[3].replaceAll(",", ""));
					correl = toFloat(data[4]);
					fb.addTC(name, tr, sd, nohsnp, correl);

				}else if (data != null && (data[0].equals("tc used"))) {

					fb.settr(Float.parseFloat(data[1]));
					fb.setTcflg(Integer.parseInt(data[2]));
					
				}else if (data != null && (data[0].equals("ploidy"))) {
					fb.getTCBean().ploidy = Float.parseFloat(data[1]);
				}else if (data != null && data[0].startsWith("chr")) {

					CNAInterval cnai = new CNAInterval(data);
					l.add(cnai);
				}else if (data != null && !data[0].startsWith("n=")
						&& !data[0].startsWith("chr")
						&& !data[0].startsWith("#")&& !data[0].startsWith("somatic")) {

					String key = data[0];					
					readstats.put(key, data);
					
				}

			}
			fb.cnalist = l;
			fb.datamap = readstats;

		}
	}

	private int toInt(String s) {

		try {
			return Integer.parseInt(s);
		} catch (Exception ex) {
		}
		return 0;
	}

	private float toFloat(String s) {

		try {
			return Float.parseFloat(s);
		} catch (Exception ex) {
		}
		return 0f;
	}

	public Filebean getBean(String sid) {

		for (Filebean fb : list) {
			if (fb.id.equals(sid)) {
				return fb;
			}
		}
		return null;
	}

}
