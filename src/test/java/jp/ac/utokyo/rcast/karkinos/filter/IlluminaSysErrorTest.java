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

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import java.io.File;
import java.io.IOException;

public class IlluminaSysErrorTest {
  public static void main(String[] arg) throws IOException {
    String twobitref = "/GLUSTER_DIST/data/Genomes/hg19_all/hg19.2bit";
    // String str
    // ="ATAAATTTAAAAGGGAAACTAATTTGGAAATCAGAAAACCACTAAGGAATTTGGGAATTAGGCTTCTGCTGCCCTCTCTGC";
    // boolean b = checkTypicalError('A','C',str);
    TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
    String chrom = "chr1";
    int pos = 155287631;
    String before10 = tgr.getGenomicSeq(chrom, pos - 11, pos - 1, true);
    String after10 = tgr.getGenomicSeq(chrom, pos + 1, pos + 11, true);
    System.out.println(before10 + "\t" + after10);
    String str = tgr.getGenomicSeq(chrom, pos - 40, pos + 40, true);
    String before2 = IlluminaSysError.getBf(before10, 2);
    String after2 = IlluminaSysError.getAfter(after10, 2);
    System.out.println(before2 + "\t" + after2);
    int polynuc = IlluminaSysError.polynuc('G', before10, after10);
    System.out.println(polynuc);
    boolean pa1 = IlluminaSysError.polyAorT(before10, 6);
    boolean pa2 = IlluminaSysError.polyAorT(after10, 6);
    if (pa1 || pa2) {
      System.out.println("true");
    }
    // boolean b = checkTypicalError('A','C',str);
  }
}
