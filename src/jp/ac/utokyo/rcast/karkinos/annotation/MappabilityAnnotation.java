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
package jp.ac.utokyo.rcast.karkinos.annotation;

import java.io.IOException;

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

public class MappabilityAnnotation {
	
	String file = "";
	BBFileReader bbReader = null;
	boolean init = false;
	public MappabilityAnnotation(String s) throws IOException{
		if(s!=null){
		 file = s;
		 bbReader = new BBFileReader(s);
		 init = true;
		}
	}
	
	public boolean isInit() {
		return init;
	}

	public float getMappability(String chrom,int pos){
		
		float ret = -1;
		BigWigIterator bi = bbReader.getBigWigIterator(chrom, pos, chrom, pos, false);
		if(bi!=null&&bi.hasNext()){
			WigItem wig = bi.next();
			ret = wig.getWigValue();
		}
		return ret;
	}
	

}
