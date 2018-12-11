package jp.ac.utokyo.rcast.karkinos.exec;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.Test;
import java.io.File;
import java.util.SortedMap;
import java.util.TreeMap;
import jp.ac.utokyo.rcast.karkinos.utils.Interval;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

// TODO: Add more tests.
public class CaptureHolderTest {
  private static final File bed = new File("test-resources/bed/test.bed");
  private static final File twobit = new File("test-resources/twobit/test.2bit");

  private static boolean isSameCapInterval(final CapInterval a, final CapInterval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && a.isGene() == b.isGene();
  }

  @Test
  void loadTargetBedTest() {
    final CaptureHolder ch = new CaptureHolder();
    try {
      ch.loadTargetBedFirstTime(bed, new TwoBitGenomeReader(twobit));
    } catch (final Exception e) {
      e.printStackTrace();
    }

    //        12345678901234567890123456789012
    // 1,11   ===========
    // 21,31                      ===========
    assertFalse(ch.getMap().isEmpty());
    assertTrue(ch.getMap().containsKey("chr1"));
    final TreeMap<Integer, CapInterval> m = ch.getMap().get("chr1");
    assertEquals(2, m.size());
    assertTrue(isSameCapInterval(new CapInterval("chr1", 1, 11, false), m.get(1)));
    assertTrue(isSameCapInterval(new CapInterval("chr1", 21, 31, false), m.get(21)));
  }

  @Test
  void getIntersectCapintervalTest() {
    final CaptureHolder ch = new CaptureHolder();
    try {
      ch.loadTargetBedFirstTime(bed, new TwoBitGenomeReader(twobit));
    } catch (final Exception e) {
      e.printStackTrace();
    }

    //        12345678901234567890123456789012
    // 1,11   ===========
    // 21,31                      ===========
    // 1,11   -----------
    // 1,20   --------------------
    // 1,21   ---------------------
    // 12,20             ---------
    {
      final SortedMap<Integer, CapInterval> m = ch.getIntersectCapinterval(new Interval("chr1", 1, 11));
      assertEquals(1, m.size());
      assertTrue(isSameCapInterval(new CapInterval("chr1", 1, 11, false), m.get(1)));
    }

    {
      final SortedMap<Integer, CapInterval> m = ch.getIntersectCapinterval(new Interval("chr1", 1, 20));
      assertEquals(1, m.size());
      assertTrue(isSameCapInterval(new CapInterval("chr1", 1, 11, false), m.get(1)));
    }

    {
      final SortedMap<Integer, CapInterval> m = ch.getIntersectCapinterval(new Interval("chr1", 1, 21));
      assertEquals(2, m.size());
      assertTrue(isSameCapInterval(new CapInterval("chr1", 1, 11, false), m.get(1)));
      assertTrue(isSameCapInterval(new CapInterval("chr1", 21, 31, false), m.get(21)));
    }

    assertTrue(ch.getIntersectCapinterval(new Interval("chr1", 12, 20)).isEmpty());
  }
}
