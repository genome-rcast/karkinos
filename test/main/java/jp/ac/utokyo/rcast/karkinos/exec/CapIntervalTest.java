package jp.ac.utokyo.rcast.karkinos.exec;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.util.CloseableIterator;
import java.io.File;

public class CapIntervalTest {
  private static final double delta = 1e-4;

  private static boolean isSameCapInterval(final CapInterval a, final CapInterval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && Math.abs(a.getCgParcent() - b.getCgParcent()) < delta
      && Math.abs(a.getDuality() - b.getDuality()) < delta;
  }

  @Test
  void getDepthTest() {
    final CapInterval ci = new CapInterval("chr1", 1, 10, false);
    assertEquals(ci.getDepth(20), 2f, delta);
  }

  @Test
  void mergeTest() {
    CapInterval ci = new CapInterval("chr1", 2, 11, false);

    ci.merge(new CapInterval("chr1", 3, 12, false));
    assertTrue(isSameCapInterval(new CapInterval("chr1", 2, 12, false, 0.0f, 1.8181f), ci));

    ci.merge(new CapInterval("chr1", 1, 10, false));
    assertTrue(isSameCapInterval(new CapInterval("chr1", 1, 12, false, 0.0f, 2.5f), ci));
  }

  @Test
  void intersectTest() {
    final CapInterval ci = new CapInterval("chr1", 5, 10, false);

    //      12345678901234567890
    // ci       ======
    // arg      ------
    assertTrue(ci.intersect(ci));

    //      12345678901234567890
    // ci       ======
    // arg  ----
    assertFalse(ci.intersect(new CapInterval("chr1", 1, 4, false)));

    //      12345678901234567890
    // ci       ======
    // arg  -----
    assertTrue(ci.intersect(new CapInterval("chr1", 1, 5, false)));

    //      12345678901234567890
    // ci       ======
    // arg           -----
    assertTrue(ci.intersect(new CapInterval("chr1", 10, 14, false)));

    //      12345678901234567890
    // ci       ======
    // arg            ----
    assertFalse(ci.intersect(new CapInterval("chr1", 11, 14, false)));
  }

  @Test
  void intersectFromSAMFileTest() {
    //           123456789012345678901234567890123456789012345
    // ci 10,15           ======
    //    7,22         ----------------
    //    9,18           ----------
    //    9,14           ------
    //    16,40                 -------------------------
    //    29,33                              -----
    //    37,45                                      ---------
    final File file = new File("test-resources/sam/test.sam");

    final CapInterval ci = new CapInterval("chr1", 10, 15, false);
    final SAMFileReader rdr = new SAMFileReader(file);
    final CloseableIterator<SAMRecord> it = rdr.iterator();
    final boolean[] expected = new boolean[]{true, true, true, false, false, false};
    int i = 0;
    while (it.hasNext()) {
      final SAMRecord record = it.next();
      assertEquals(expected[i], ci.intersect(record));
      ++i;
    }
    it.close();
  }

  @Test
  void includeTest() {
    final CapInterval ci = new CapInterval("chr1", 5, 10, false);

    //      12345678901234567890
    // ci       ======
    // pos     .
    assertFalse(ci.include(4));

    //      12345678901234567890
    // ci       ======
    // pos      .
    assertTrue(ci.include(5));

    //      12345678901234567890
    // ci       ======
    // pos           .
    assertTrue(ci.include(10));

    //      12345678901234567890
    // ci       ======
    // pos            .
    assertFalse(ci.include(11));
  }
}
