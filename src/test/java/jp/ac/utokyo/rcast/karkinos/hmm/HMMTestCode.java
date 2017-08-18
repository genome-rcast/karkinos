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
package jp.ac.utokyo.rcast.karkinos.hmm;

import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixture;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussianMixtureFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfInteger;
import be.ac.ulg.montefiore.run.jahmm.OpdfIntegerFactory;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;

public class HMMTestCode {
	
public static void main(String[] aa){
		
	
	                                                 

	OpdfGaussianFactory factory = new OpdfGaussianFactory();
	Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(4, factory);
     
	hmm.setPi(0, 0.05);
	hmm.setPi(1, 0.90);
	hmm.setPi(2, 0.04);
	hmm.setPi(3, 0.01);


	hmm.setOpdf(0, new OpdfGaussian(0.5,0.11));
	hmm.setOpdf(1, new OpdfGaussian(1,0.11));
	hmm.setOpdf(2, new OpdfGaussian(1.5,0.11));
	hmm.setOpdf(3, new OpdfGaussian(2,0.1));
		

		hmm.setAij(0, 1, 0.1);
		hmm.setAij(0, 0, 0.9);
		hmm.setAij(0, 2, 0);
		hmm.setAij(0, 3, 0);
		
		hmm.setAij(1, 0, 0.05);
		hmm.setAij(1, 1, 0.9);
		hmm.setAij(1, 2, 0.05);
		hmm.setAij(1, 3, 0);
		
		hmm.setAij(2, 0, 0);
		hmm.setAij(2, 1, 0.08);
		hmm.setAij(2, 2, 0.9);
		hmm.setAij(2, 3, 0.2);
		
		hmm.setAij(3, 0, 0);
		hmm.setAij(3, 1, 0);
		hmm.setAij(3, 2, 0.1);
		hmm.setAij(3, 3, 0.9);
		
		List<ObservationReal>  l = new ArrayList<ObservationReal>();
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(0.5));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(0));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(1));
		l.add(new ObservationReal(2));
		l.add(new ObservationReal(2));
		l.add(new ObservationReal(2));
		l.add(new ObservationReal(2));
		l.add(new ObservationReal(2));
		l.add(new ObservationReal(2));

		int[] ary = hmm.mostLikelyStateSequence(l);
		for(int n:ary){
			System.out.println(n);
		}
		
	}

}
