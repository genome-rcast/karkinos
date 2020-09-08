package jp.ac.utokyo.rcast.karkinos.exec;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

public class CopyNumberIntervalTest {
  @Test
  void lengthTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {{
      setStart(2);
      setEnd(4);
    }};
    assertEquals(3, cni.length());
  }
}
