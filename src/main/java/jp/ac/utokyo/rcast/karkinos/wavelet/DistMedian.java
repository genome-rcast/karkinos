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
package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

public class DistMedian implements java.io.Serializable{

	public float MIN = 0.25f;
	public float MAX = 4.5f;
	public static final float interval = 0.025f;
	public boolean considermean = false;

	Map<Integer, SummaryStatistics> map = new LinkedHashMap<Integer, SummaryStatistics>();
	SummaryStatistics nettotal = new SummaryStatistics();
	
	public DistMedian(double min, double max,boolean considermean) {
		MIN = (float)min;
		MAX = (float)max;
		this.considermean = considermean;
	}

	public DistMedian() {
		
	}

	public void addValue(double z) {

		if ((z < MIN) || (z > MAX) || (z == 1)) {
			return;
		}
		nettotal.addValue(z);
		
		int n = getIntervalIndex(z);
		SummaryStatistics ss = null;
		if (map.containsKey(n)) {
			ss = map.get(n);
		} else {
			ss = new SummaryStatistics();
			map.put(n, ss);
		}
		ss.addValue(z);
	}

	private int getIntervalIndex(double z) {

		//
		int n = (int) (z * 1000) / (int) (interval * 1000);
		return n;

	}

	public double getDistributionMedian() {
		
		double defval = nettotal.getGeometricMean();
		try{
			double dm = _getDistributionMedian();
			if(considermean&&Math.abs(defval-dm)>0.1){
				return dm;
			}
			
		}catch(Exception ex){
			
		}
		//return default val
		return defval;
	}
	
	public double _getDistributionMedian() {

		// get max distributed value
		Iterator<Integer> ite = map.keySet().iterator();
		//
		int max = 0;
		int maxidx = 0;
		int secondbestidx = 0;

		// get Max distributed point
		while (ite.hasNext()) {

			int idx = ite.next();
			SummaryStatistics ss = map.get(idx);
			//
			if (ss.getN() > max) {
				max = (int) ss.getN();
				secondbestidx = maxidx;
				maxidx = idx;
			}

		}
		int neighbor = 20;

		List<Point2D.Double> datapoints = new ArrayList<Point2D.Double>();
		maxidx = getNearToMean(maxidx,secondbestidx);

		// get +-5 max point
		for (int n = maxidx - neighbor; n <= maxidx + neighbor; n++) {

			SummaryStatistics ss = map.get(n);
			if (ss != null) {

				//
				datapoints.add(new Point2D.Double(ss.getMean(), ss.getN()));

			}

		}

		int n = 0;
		double X = 0;
		double X2 = 0;
		double X3 = 0;
		double X4 = 0;

		double Y = 0;
		double XY = 0;
		double X2Y = 0;

		// LSM to find 2nd degree eq
		for (Point2D.Double point : datapoints) {

			//
			double x0 = point.getX();
			double y0 = point.getY();
			double x2 = Math.pow(x0, 2);
			double x3 = Math.pow(x0, 3);
			double x4 = x2 * x2;

			X += x0;
			X2 += x2;
			X3 += x3;
			X4 += x4;
			n++;
			
			Y += y0;
			XY += x0 * y0;
			X2Y += x2 * y0;

		}
		// resolve equation
		double[][] vals = { { n, X, X2 }, { X, X2, X3 }, { X2, X3, X4 } };
		double[] rhs = { Y, XY, X2Y };

		RealMatrix a = new Array2DRowRealMatrix(vals);

		double a1 = 0;
		double b1 = 0;
		double c1 = 0;

		try {
			DecompositionSolver solver = new LUDecompositionImpl(a).getSolver();
			RealVector b = new ArrayRealVector(rhs);
			RealVector x = solver.solve(b);
			a1 = x.getEntry(2);
			b1 = x.getEntry(1);
			c1 = x.getEntry(0);
		} catch (Exception ex) {
			//ex.printStackTrace();
		}

		if (a1 >= 0) {
			// oops, something is wrong
			//return max interval mean for next to best efort
			if(map.containsKey(maxidx)){
				return map.get(maxidx).getMean();
			}else{
				return nettotal.getGeometricMean();
			}
			
		}

		// derivative to find max
		double maxX = -1 * (b1 / (2 * a1));

		//

		return maxX;
	}

	private int getNearToMean(int maxidx, int secondbestidx) {
		
		if(map.containsKey(secondbestidx)){
			
			double netmean = nettotal.getMean();
			double mean1 = map.get(maxidx).getMean();
			double mean2 = map.get(secondbestidx).getMean();
			boolean oneIsclose = Math.abs(netmean-mean1)<Math.abs(netmean-mean2);
			return oneIsclose?maxidx:secondbestidx;
		}
		return maxidx;
	}

}
