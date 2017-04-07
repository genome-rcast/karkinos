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
package jp.ac.utokyo.rcast.karkinos.readssummary;

public class GeneExonsTest {
  public static void main(String[] arg) {
    String ref = "/GLUSTER_DIST/data/Genomes/ucscgene_hg19/totoki_hg19/refFlat_120402.txt";
    GeneExons ex = new GeneExons(ref);
  }
}
