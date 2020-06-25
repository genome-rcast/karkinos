package jp.ac.utokyo.rcast.karkinos.readssummary;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

public class IntervalTest {
  private static boolean isSameInterval(final Interval a, final Interval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && a.getDepth() == b.getDepth()
      && a.getRefseqid().equals(b.getRefseqid())
      && a.getGeneSymbol().equals(b.getGeneSymbol());
  }

  @Test
  void extendIntervalTest() {
    final Interval iv = new Interval("chr1", 2, 4, "NR_037415", "MIR3620");
    iv.depth = 0;

    //      12345
    // iv    ===
    // pos  .
    //  =>   ===  doesn't extend
    assertFalse(iv.extendInterval("chr1", 1, 0));

    //      12345
    // iv    ===
    // pos     .
    //  =>   ===  doesn't extend
    assertFalse(iv.extendInterval("chr1", 4, 0));

    //      12345
    // iv    ===
    // pos      .
    //  =>   ====  extend
    assertTrue(iv.extendInterval("chr1", 5, 0));
    final Interval expected = new Interval("chr1", 2, 5, "NR_037415", "MIR3620");
    expected.depth = 0;
    assertTrue(isSameInterval(expected, iv));
  }

  private static Stream<Arguments> containArgs() {
    return Stream.of(
        Arguments.of(false, "chr1", 1),
        Arguments.of(true, "chr1", 2),
        Arguments.of(true, "chr1", 4),
        Arguments.of(false, "chr1", 5),
        Arguments.of(false, "chr2", 2));
  }

  @ParameterizedTest
  @MethodSource("containArgs")
  void containTest(final boolean expected, final String chr, final int pos) {
    final Interval iv = new Interval("chr1", 2, 4, "NR_037415", "MIR3620");
    assertEquals(expected, iv.contain(chr, pos));
  }
}
