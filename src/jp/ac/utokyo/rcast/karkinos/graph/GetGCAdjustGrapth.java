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
package jp.ac.utokyo.rcast.karkinos.graph;

import java.awt.image.BufferedImage;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.tc33.jheatchart.HeatChart;


import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.FunctionRegression;

public class GetGCAdjustGrapth {

	private final static int SIZE = 2;
	public static List<DisplayObject> getChartLists(DataSet dataset) {

		FunctionRegression fr = dataset.getFunctionRegression();
		List<DisplayObject> list = new ArrayList<DisplayObject>();
		//
		List oList = new ArrayList();
		oList.add(getMeanGraph(fr));
		oList.add(getSDGraph(fr));
		list.add(new DisplayObject(oList, SIZE, "GC% Adjustment"));
		return list;
	}
	
	private static Object getSDGraph(FunctionRegression fr) {
		
		double[][] sd = fr.getSDAry();
		HeatChart map = new HeatChart(sd); 
		
		map.setTitle("s.d.("+format(map.getLowValue())+" to "+format(map.getHighValue())+")");
		map.setYAxisLabel("Bait length");
		map.setXAxisLabel("GC %");
		map.setXValues(fr.getCGLabel());
		//map.setXValuesHorizontal(true);
		map.setYValues(fr.getBaitLabel());
		java.awt.Image im = map.getChartImage();
		return im;
	}

	private static java.awt.Image getMeanGraph(FunctionRegression fr) {
		
		
		double[][] mean = fr.getMeanAry();
		HeatChart map = new HeatChart(mean); 
		map.setTitle("Mean T/N depth ratio ("+format(map.getLowValue())+" to "+format(map.getHighValue())+")");
		map.setYAxisLabel("Bait length");
		map.setXAxisLabel("GC %");	
		map.setXValues(fr.getCGLabel());
		map.setYValues(fr.getBaitLabel());
		//map.setXValuesHorizontal(true);
		java.awt.Image im = map.getChartImage();
		return im;
		
	}
	
	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}
	
}
