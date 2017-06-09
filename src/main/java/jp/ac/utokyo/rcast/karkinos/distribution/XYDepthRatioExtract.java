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
package jp.ac.utokyo.rcast.karkinos.distribution;

import java.text.NumberFormat;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.utils.CalcUtils;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class XYDepthRatioExtract {

	DataSet dataset;

	public XYDepthRatioExtract(DataSet dataset) {
		this.dataset = dataset;
	}

	public XYSeriesCollection getTumorRatioDepth(int i) {

		XYSeriesCollection data = new XYSeriesCollection();
		Map<Float, XYSeries> map = new LinkedHashMap<Float, XYSeries>();
		for (float degree = 0.5f; degree <= 2; degree = degree + 5f) {
			//
			String ns = format((degree * 2));
			XYSeries n = new XYSeries("n=" + ns);
			map.put(degree, n);

		}
		for (float degree = 0.5f; degree <= 2; degree = degree + 0.25f) {
			//
			String ns = format((degree * 2));
			if (!map.containsKey(degree)) {
				XYSeries n = new XYSeries("n=" + ns);
				map.put(degree, n);
			}

		}

		for (SNVHolder snv : dataset.getSnvlist()) {

			float f = (float) snv.getCi().getVaridateVal();
			float ratio = snv.getTumor().getRatio();
			int depth = snv.getTumor().getTotalcnt();
			XYSeries xys = map.get(f);
			if(xys==null)continue;

			//
			if (i == 0) {
				if (snv.isHetroSNP()) {
					xys.add(ratio, depth);
				}
			} else if (i == 1) {
				if (snv.getFlg() == PileUP.SomaticMutation) {
					xys.add(ratio, depth);
				}
			} else if (i == 2) {
				if (snv.getFlg() == PileUP.SomaticMutation) {
					if (snv.getFilterResult() != null
							&& snv.getFilterResult().isPassFilter()) {
						xys.add(ratio, depth);
					}
				}

			} else if (i == 3) {
				if (snv.getFlg() == PileUP.SomaticMutation) {
					if (snv.getFilterResult() != null
							&& snv.getFilterResult().isPassFilter()) {
						if (CalcUtils.pass2(snv.getFilterResult())) {
							xys.add(ratio, depth);
						}
					}
				}

			}

		}

		Iterator<Float> ite = map.keySet().iterator();
		while (ite.hasNext()) {

			//
			XYSeries xys = map.get(ite.next());
			if (xys.getItemCount() > 0) {
				data.addSeries(xys);
			}
		}
		return data;
	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

}
