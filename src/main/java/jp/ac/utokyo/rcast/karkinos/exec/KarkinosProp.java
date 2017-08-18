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
package jp.ac.utokyo.rcast.karkinos.exec;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.util.Properties;

public class KarkinosProp implements java.io.Serializable {

	
	public static final float TNQdiff = 10.0f;
	
	static Properties config = null;

	public static void load(String pathToConfig) {

		try {

			InputStream is = new FileInputStream(new File(pathToConfig));
			config = new Properties();
			config.load(is);
			loadEach();

		} catch (Exception e) {
			System.err.println("fail to read property file");
			e.printStackTrace();
		}
	}

	public static int BINBITSIZE = 18;
	public static double denozeToSD = 0.02;
	public static int maxdenoiseLevel = 8;
	//min reads counts for target region when calculating
	//normal-tumor depth ratio for target region for CNV analysis
	public static int mincoverbase = 1000;
	public static int pairReadsMaxdist = 1000;

	public static int mindepth = 8;
	public static int mindepthNormal = 5;
	
	public static float maxnormalratio = 0.025f;
	public static float min_initial_tumorratio = 0.05f;
	
	public static float mintumorratioForFilter1 = 0.1f;
	public static float mintumorratio = 0.08f;
	public static float mintumorratioForResqued = 0.2f;
	public static float normalSNPthres = 0.05f;
	public static int minsupportreads = 3;

	public static float hetroSNPMin = 0.38f;
	public static float hetroSNPMax = 0.62f;

	public static float minPhredQual = 60f;
	public static float minMappability = 0.3f;
	public static float minEntropy = 0.8f;
	public static float minMisMatchRate = 0.025f;
	// deplicated not using
	public static int sample_max_bait_size = 200;

	public static double Fisher_Thres_For_Reads_Direction = 0.03;
	public static double Fisher_Thres_For_SNV_Detection = 0.15;
	
	public static final double Fisher_Thres_For_Reads_Direction2 = 0.2;
	public static final double Fisher_Thres_For_Reads_Direction3 = 0.4;
	
	
	public static float minMapQ = 20;
	public static float pvalforNoisepeak = 0.1f;
	
	public static float minPhredQualForEach = 10f;
	public static double LogtThres = -0.5;
	public static double LognThres = 2.0;
	public static int nearindelbt = 5;	
	public static int baitmergin = 300;
	public static int baysianFilterdepth = 40;
	public static int tcFilterdepth = 100;
	public static int entropyDepth = 50;
	public static float lowQualRatiothres =  0.1f;
	
	public static int unittumorread = 20;
	public static float tcFilterfreq = 0.35f;
	public static int low_normal_depth_thresLow = 20;
	public static int low_normal_depth_thresHigh = 40;
	//public static float mintumorratioOrg = 0.045f;
	public static double minTumorNormalRatio = 0.15;
	public static double minSecondMutationRelativeRatio = 0.2;	
	public static double mintumorratioAdj = 0.15;
	
	
	public static float mintumorpurity = 0.12f;
	
	public static float falseReadMismatchratio = 0.03f;
	public static float falseReadratio = 0.3f;
	public static float falseReadratio2 = 0.2f;
	
	
	public static String KEY_BINBITSIZE = "BINBITSIZE";
	public static String KEY_denozeToSD = "denozeToSD";
	public static String KEY_minsupport = "minsupportreads";
	public static String KEY_maxdenoiseLevel = "maxdenoiseLevel";
	public static String KEY_mincoverbase = "mincoverbase";
	public static String KEY_minPhredQualForEach = "minPhredQualForEach";

	public static String KEY_mindepth = "mindepth";
	public static String KEY_mindepth_normal = "mindepth_normal";
	public static String KEY_mintumorratio = "mintumorratio";
	public static String KEY_min_initial_tumorratio = "min_initial_tumorratio";
	public static String KEY_maxnormalratio = "maxnormalratio";
	public static String KEY_normalSNPthres = "normalSNPthres";

	public static String KEY_hetroSNPMin = "hetroSNPMin";
	public static String KEY_hetroSNPMax = "hetroSNPMax";
	public static String KEY_minPhredQual = "minPhredQual";
	public static String KEY_minMappability = "minMappability";

	public static String KEY_minEntropy = "minEntropy";
	public static String KEY_minMisMatchRate = "minMisMatchRate";
	public static String KEY_Fisher_Thres_For_Reads_Direction = "Fisher_Thres_For_Reads_Direction";
	public static String KEY_Fisher_Thres_For_SNV_Detection = "Fisher_Thres_For_SNV_Detection";

	public static String KEY_LognThres  ="LognThres";
	public static String KEY_LogtThres  ="LogtThres";
	public static String KEY_nearindelbt  ="nearindelbt";
	public static String KEY_baysianFilterdepth  ="baysianFilterdepth";

	


	private static void loadEach() {

		//
		BINBITSIZE = getIntProperty(KEY_BINBITSIZE,BINBITSIZE);
		denozeToSD = getFloatProperty(KEY_denozeToSD,(float)denozeToSD);
		maxdenoiseLevel = getIntProperty(KEY_maxdenoiseLevel,maxdenoiseLevel);
		mincoverbase = getIntProperty(KEY_mincoverbase,mincoverbase);

		minsupportreads = getIntProperty(KEY_minsupport,minsupportreads);
		
		mindepth = getIntProperty(KEY_mindepth,mindepth);
		mindepthNormal = getIntProperty(KEY_mindepth_normal,mindepthNormal);
		
		mintumorratio = getFloatProperty(KEY_mintumorratio,mintumorratio);
		min_initial_tumorratio
			= getFloatProperty(KEY_min_initial_tumorratio,min_initial_tumorratio);
		
		maxnormalratio = getFloatProperty(KEY_maxnormalratio,maxnormalratio);
		normalSNPthres = getFloatProperty(KEY_normalSNPthres,normalSNPthres);

		hetroSNPMin = getFloatProperty(KEY_hetroSNPMin,hetroSNPMin);
		hetroSNPMax = getFloatProperty(KEY_hetroSNPMax,hetroSNPMax);

		minPhredQual = getFloatProperty(KEY_minPhredQual,minPhredQual);
		minMappability = getFloatProperty(KEY_minMappability,minMappability);
		minEntropy = getFloatProperty(KEY_minEntropy,minEntropy);
		//minMisMatchRate = getFloatProperty(KEY_minMisMatchRate,minMisMatchRate);
		
		Fisher_Thres_For_Reads_Direction = getDoubleProperty(KEY_Fisher_Thres_For_Reads_Direction,Fisher_Thres_For_Reads_Direction);
		Fisher_Thres_For_SNV_Detection = getDoubleProperty(KEY_Fisher_Thres_For_SNV_Detection,Fisher_Thres_For_SNV_Detection);
		minPhredQualForEach =  getFloatProperty(KEY_minPhredQualForEach,minPhredQualForEach);
		
		LognThres = getFloatProperty(KEY_LognThres,(float)LognThres);
		//LogtThres = getFloatProperty(KEY_LogtThres,(float)LogtThres);
		nearindelbt = getIntProperty(KEY_nearindelbt,nearindelbt);
		baysianFilterdepth = getIntProperty(KEY_baysianFilterdepth,baysianFilterdepth);

		
	}
	
	public static String getInfoString(){
		
		StringBuffer sb = new StringBuffer();
		KarkinosProp inst = new KarkinosProp();
		for(Field f:inst.getClass().getFields()){
			try {
				sb.append(f.getName()+"="+f.get(inst));
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				
			} catch (IllegalAccessException e) {
				// TODO Auto-generated catch block
			}
		}
		return sb.toString();
		
	}

	public static int getIntProperty(String key,int n) throws NumberFormatException {

		try{
			return Integer.parseInt(config.getProperty(key));
		}catch(Exception ex){
			
			System.err.println("could not find key="+key +" in property file \n" +
					"default val "+n+" is used");
			return n;
		}
	}

	public static double getDoubleProperty(String key,double d) {

		try{
			return Double.parseDouble(config.getProperty(key));
		}catch(Exception ex){
			System.err.println("coud not find key="+key +" in property file \n" +
					"default val "+d+" is used");
			return d;
		}
	}
	

	public static float getFloatProperty(String key,float f) {
		
		try{
			return Float.parseFloat(config.getProperty(key));
		}catch(Exception ex){
			System.err.println("coud not find key="+key +" in property file \n" +
					"default val "+f+" is used");
			return f;
		}
	}

}
