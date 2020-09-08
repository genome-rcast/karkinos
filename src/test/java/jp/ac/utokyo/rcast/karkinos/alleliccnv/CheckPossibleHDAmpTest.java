package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;

// TODO: This test class only aims at testing CopyNumberInterval and CapInterval.
// We should add more tests in the future.
public class CheckPossibleHDAmpTest {
  private static final int unit = 1000000;

  @Test
  void in1mTest() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
    final Method in1m = CheckPossibleHDAmp.class.getDeclaredMethod("in1m",
        CopyNumberInterval.class, CapInterval.class);
    in1m.setAccessible(true);

    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {{
      setStart(2 * unit + 4);
      setEnd(2 * unit + 6);
    }};

    //      1234567890
    // cni     ===
    // arg ----
    assertFalse((boolean)in1m.invoke(CheckPossibleHDAmp.class, cni,
          new CapInterval("chr1", unit + 1, unit + 3, false)));

    //      123456789
    // cni     ===
    // arg -----
    assertTrue((boolean)in1m.invoke(CheckPossibleHDAmp.class, cni,
          new CapInterval("chr1", unit + 1, unit + 4, false)));

    //      123456789
    // cni     ===
    // arg       -----
    assertTrue((boolean)in1m.invoke(CheckPossibleHDAmp.class, cni,
          new CapInterval("chr1", 3 * unit + 6, 3 * unit + 7, false)));

    //      123456789
    // cni     ===
    // arg        ----
    assertFalse((boolean)in1m.invoke(CheckPossibleHDAmp.class, cni,
          new CapInterval("chr1", 3 * unit + 7, 3 * unit + 8, false)));
  }
}
