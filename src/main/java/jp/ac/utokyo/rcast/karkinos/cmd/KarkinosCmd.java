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
package jp.ac.utokyo.rcast.karkinos.cmd;

import jp.ac.utokyo.rcast.karkinos.exec.TumorGenotyper;

import java.util.Arrays;
import java.util.Optional;

public class KarkinosCmd {
  public static void main(String[] arg) throws Exception {
    if (arg == null || arg.length == 0 || arg[0] == null) {
      printMessage();
      return;
    }

    final String[] arg2 = Arrays.copyOfRange(arg,1,arg.length);
    switch (arg[0]) {
      case "analysis":
        TumorGenotyper.main(arg2);
        break;
      default:
        printMessage();
        break;
    }
  }

  public static Optional<String> getVersion() {
    final String version = KarkinosCmd.class.getPackage().getImplementationVersion();
    final String title = KarkinosCmd.class.getPackage().getImplementationTitle();
    if (version != null && title != null) {
      return Optional.of(title + " " + version);
    }
    return Optional.empty();
  }

  private static void printMessage() {
    getVersion().ifPresent(System.err::println);
    System.err.println("usage: karkinos.jar <command> options");
    System.err.println("Command: ");
    System.err.println("  analysis	analysis SNV,CSV,tumor rate from bamfiles");
  }
}
