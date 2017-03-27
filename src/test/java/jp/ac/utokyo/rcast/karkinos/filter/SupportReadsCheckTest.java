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
package jp.ac.utokyo.rcast.karkinos.filter;

import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.exec.IndelInfo;
import jp.ac.utokyo.rcast.karkinos.exec.PileUP;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class SupportReadsCheckTest {
  public static void main(String[] arg) {
    //
    // String bam =
    // "/GLUSTER_DIST/data/users/yamamoto/exome/OvCa/karkinos4.1.8/summary_cfDNA1/bam/Se-67-tumor25-Se-67_b-DNA_tumor_genome.bam";
    // String nbam =
    // "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT30-2-MT30N/normal/MT30-2-MT30N_normal_genome.bam";
    // String bam =
    // "/GLUSTER_DIST/data/users/yamamoto/exome/Glioma/MT/MT30-2-MT30N/tumor/MT30-2-MT30N_tumor_genome.bam";
    String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
    String nbam = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_normal_genome.bam";
    String bam = "/data3/users/yamamoto/exome/CRC/karkinos4.1.11/summary_CRC_all_samples/bam/CRC107_T-CRC107_N_tumor_genome.bam";
    String middelfile = "/GLUSTER_DIST/data/users/yamamoto/exome/CRC/karkinos2.0.3/CRC_107_T-CRC_107_N/sobj";
    TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
    DbSNPAnnotation dbAnno = null;
    SupportReadsCheck inst = new SupportReadsCheck(nbam, bam, tgr, dbAnno);
    //
    PileUPResult pileUPResult = new PileUPResult();
    pileUPResult.setGenomeRef('C');
    IndelInfo ii = new IndelInfo();
    // ii.indel = true;
    // ii.length = 17;
    // ii.cnt = 13;
    // // ii.insersion = "G";
    pileUPResult.setBaseAndQual('A', (byte) 50, 1, ii);
    Map<String, Integer> snppos = new HashMap<String, Integer>();
    try {
      inst.checkSupportReads("chr12", 25398284, pileUPResult, PileUP.SomaticMutation, 0.2f,
                             PileUP.SomaticMutation, false, 10, snppos, 0.2f, 0.2f, 0f, 0f);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}
