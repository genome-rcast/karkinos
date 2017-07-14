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

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

public class CaptureHolderTest {
  public static void main(String[] arg) {
    CaptureHolder inst = new CaptureHolder();
    //String bed = "/GLUSTER_DIST/data/users/ueda/testframework/script/vcrome2.1.bed";
    String bed = "/GLUSTER_DIST/data/Genomes/karkinos/genome/SureSelectV5plusLincRNA_Regions.bed";
    String tb = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
    TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(tb));
    try {
      inst.loadTargetBed(bed, tgr);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}
