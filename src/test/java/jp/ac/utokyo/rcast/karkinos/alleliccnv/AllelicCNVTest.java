package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;

// TODO: This test class only aims at testing CopyNumberInterval and CapInterval.
// We should add more tests in the future.
public class AllelicCNVTest {
  private static final double delta = 1e-4;

  private static boolean isSameCopyNumberInterval(final CopyNumberInterval a, final CopyNumberInterval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && a.getAaf() == b.getAaf()
      && a.getBaf() == b.getBaf()
      && Math.abs(a.getCopynumber() - b.getCopynumber()) < delta
      && a.isAllelic() == b.isAllelic()
      && a.isHdeletion() == b.isHdeletion()
      && a.getNoSNP() == b.getNoSNP();
  }

  private static boolean isSameList(final List<CopyNumberInterval> a, final List<CopyNumberInterval> b) {
    if (a.size() != b.size()) {
      return false;
    }
    for (int i = 0; i < a.size(); ++i) {
      if (!isSameCopyNumberInterval(a.get(i), b.get(i))) {
        return false;
      }
    }
    return true;
  }

  @Test
  void analyseTest() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
    final Method analyse = AllelicCNV.class.getDeclaredMethod("analyse", List.class, List.class);
    analyse.setAccessible(true);

    // 123456
    //  ===
    //   ===
    final List<CopyNumberInterval> allelicLOHLow = new ArrayList<CopyNumberInterval>() {{
      final CopyNumberInterval c1 = new CopyNumberInterval("chr1") {{
        setStart(2);
        setEnd(4);
        setCopynumber(1);
      }};
      add(c1);

      final CopyNumberInterval c2 = new CopyNumberInterval("chr1") {{
        setStart(3);
        setEnd(5);
        setCopynumber(1);
      }};
      add(c2);
    }};

    // 123456
    // --
    //  --
    //   --
    //    --
    //     --
    final List<CopyNumberInterval> allelicLOHhigh = new ArrayList<CopyNumberInterval>() {{
      final CopyNumberInterval c1 = new CopyNumberInterval("chr1") {{
        setStart(1);
        setEnd(2);
        setCopynumber(1);
      }};
      add(c1);

      final CopyNumberInterval c2 = new CopyNumberInterval("chr1") {{
        setStart(2);
        setEnd(3);
        setCopynumber(1);
      }};
      add(c2);

      final CopyNumberInterval c3 = new CopyNumberInterval("chr1") {{
        setStart(3);
        setEnd(4);
        setCopynumber(1);
      }};
      add(c3);

      final CopyNumberInterval c4 = new CopyNumberInterval("chr1") {{
        setStart(4);
        setEnd(5);
        setCopynumber(1);
      }};
      add(c4);

      final CopyNumberInterval c5 = new CopyNumberInterval("chr1") {{
        setStart(5);
        setEnd(6);
        setCopynumber(1);
      }};
      add(c5);
    }};

    final List<CopyNumberInterval> actual = (List<CopyNumberInterval>)analyse.invoke(
        AllelicCNV.class, allelicLOHLow, allelicLOHhigh);

    final List<CopyNumberInterval> expected = new ArrayList<CopyNumberInterval>() {{
      final CopyNumberInterval c1 = new CopyNumberInterval("chr1") {{
        setStart(2);
        setEnd(3);
        setCopynumber(0);
        setHdeletion(true);
      }};
      add(c1);

      final CopyNumberInterval c2 = new CopyNumberInterval("chr1") {{
        setStart(3);
        setEnd(4);
        setCopynumber(0);
        setHdeletion(true);
      }};
      add(c2);

      final CopyNumberInterval c3 = new CopyNumberInterval("chr1") {{
        setStart(4);
        setEnd(5);
        setCopynumber(0);
        setHdeletion(true);
      }};
      add(c3);
    }};

    assertTrue(isSameList(expected, actual));
  }
}
