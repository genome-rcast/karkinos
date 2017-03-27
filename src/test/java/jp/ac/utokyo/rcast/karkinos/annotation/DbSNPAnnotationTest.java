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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class DbSNPAnnotationTest {
  public static void main(String[] arg) throws IOException {
    String dbSNP = "/GLUSTER_DIST/data/users/ueda/SNVtest/hg19_ALL.sites.2012_02.txt";
    FileInputStream fis = new FileInputStream(dbSNP);
    BufferedReader br = new BufferedReader(new InputStreamReader(fis));
    try {
      int totalcnt = 0;
      int init = fis.available();
      String chr = "";
      for (; ; ) {
        int point = init - fis.available();
        String line = br.readLine();
        if (line == null)
          break;
        totalcnt++;
        String[] sa = line.split("\t");
        String _chr = sa[0];
        if (!chr.equals(_chr)) {
          System.out.println(_chr + "\t" + totalcnt + "\t" + point);
        }
        chr = _chr;
        int pos = Integer.parseInt(sa[1]);
      }
    } finally {
      br.close();
    }
  }
}
