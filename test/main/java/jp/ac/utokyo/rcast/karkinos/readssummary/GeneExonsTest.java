package jp.ac.utokyo.rcast.karkinos.readssummary;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertIterableEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.Stream;

// TODO: This test class only aims at testing Interval.
// We should add more tests in the future.
public class GeneExonsTest {
  private static final String refGene = "test-resources/refgene/refGene.txt";

  private static final Map<String, Integer> expectedCounter = new HashMap<String, Integer>() {{
    put("", 9);
  }};

  private static final Map<String, TreeMap<Integer, Interval>> expectedGeneMap = new LinkedHashMap<String, TreeMap<Integer, Interval>>() {{
    put("chr1", new TreeMap<Integer, Interval>() {{
      put(228284964, new Interval("chr1", 228284964, 228285042, "NR_037415", "2326"));
    }});
    put("chr6", new TreeMap<Integer, Interval>() {{
      put(133065009, new Interval("chr6", 133065009, 133079147, "NR_034174", "1600"));
    }});
    put("chr15", new TreeMap<Integer, Interval>() {{
      put(49447956, new Interval("chr15", 49447956, 49622002, "NM_001001556", "120"));
    }});
    put("chr16", new TreeMap<Integer, Interval>() {{
      put(16400227, new Interval("chr16", 16400227, 16400291, "NR_128713", "710"));
      put(18499555, new Interval("chr16", 18499555, 18499619, "NR_128713", "726"));
      put(15001574, new Interval("chr16", 15001574, 15001638, "NR_128712", "699"));
    }});
    put("chr17", new TreeMap<Integer, Interval>() {{
      put(1617197, new Interval("chr17", 1617197, 1617281, "NR_029494", "597"));
    }});
    put("chr17_ctg5_hap1", new TreeMap<Integer, Interval>() {{
      put(592236, new Interval("chr17_ctg5_hap1", 592236, 595386, "NR_034172", "589"));
    }});
    put("chr19", new TreeMap<Integer, Interval>() {{
      put(16466055, new Interval("chr19", 16466055, 16582823, "NR_047665", "88"));
    }});
  }};

  private static final Map<String, TreeMap<Integer, Interval>> expectedMap = new LinkedHashMap<String, TreeMap<Integer, Interval>>() {{
    put("chr1", new TreeMap<Integer, Interval>() {{
      final Interval iv = new Interval("chr1", 228285042, 228285042, "NR_037415", "2326");
      iv.exonidx = 1;
      iv.depth = 228285042;
      put(228285042, iv);
    }});
    put("chr15", new TreeMap<Integer, Interval>() {{
      final Interval iv1 = new Interval("chr15", 49448193, 49448213, "NM_001001556", "120");
      iv1.exonidx = 1;
      iv1.depth = 49448213;
      put(49448193, iv1);

      final Interval iv2 = new Interval("chr15", 49493358, 49493447, "NM_001001556", "120");
      iv2.exonidx = 2;
      iv2.depth = 49493447;
      put(49493358, iv2);

      final Interval iv3 = new Interval("chr15", 49509386, 49509510, "NM_001001556", "120");
      iv3.exonidx = 3;
      iv3.depth = 49509510;
      put(49509386, iv3);

      final Interval iv4 = new Interval("chr15", 49528047, 49528138, "NM_001001556", "120");
      iv4.exonidx = 4;
      iv4.depth = 49528138;
      put(49528047, iv4);

      final Interval iv5 = new Interval("chr15", 49531417, 49531564, "NM_001001556", "120");
      iv5.exonidx = 5;
      iv5.depth = 49531564;
      put(49531417, iv5);

      final Interval iv6 = new Interval("chr15", 49574183, 49574282, "NM_001001556", "120");
      iv6.exonidx = 6;
      iv6.depth = 49574282;
      put(49574183, iv6);

      final Interval iv7 = new Interval("chr15", 49575762, 49575915, "NM_001001556", "120");
      iv7.exonidx = 7;
      iv7.depth = 49575915;
      put(49575762, iv7);

      final Interval iv8 = new Interval("chr15", 49584523, 49584734, "NM_001001556", "120");
      iv8.exonidx = 8;
      iv8.depth = 49584734;
      put(49584523, iv8);

      final Interval iv9 = new Interval("chr15", 49611800, 49612002, "NM_001001556", "120");
      iv9.exonidx = 9;
      iv9.depth = 49612002;
      put(49611800, iv9);

      final Interval iv10 = new Interval("chr15", 49620148, 49620356, "NM_001001556", "120");
      iv10.exonidx = 10;
      iv10.depth = 49620356;
      put(49620148, iv10);
    }});
    put("chr16", new TreeMap<Integer, Interval>() {{
      final Interval iv1 = new Interval("chr16", 15001638, 15001638, "NR_128712", "699");
      iv1.exonidx = 1;
      iv1.depth = 15001638;
      put(15001638, iv1);

      final Interval iv2 = new Interval("chr16", 16400291, 16400291, "NR_128713", "710");
      iv2.exonidx = 1;
      iv2.depth = 16400291;
      put(16400291, iv2);

      final Interval iv3 = new Interval("chr16", 18499619, 18499619, "NR_128713", "726");
      iv3.exonidx = 1;
      iv3.depth = 18499619;
      put(18499619, iv3);
    }});
    put("chr17", new TreeMap<Integer, Interval>() {{
      final Interval iv = new Interval("chr17", 1617281, 1617281, "NR_029494", "597");
      iv.exonidx = 1;
      iv.depth = 1617281;
      put(1617281, iv);
    }});
  }};

  private static boolean isSameInterval(final Interval a, final Interval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && a.getRefseqid().equals(b.getRefseqid())
      && a.getExonidx() == b.getExonidx()
      && a.getGeneSymbol().equals(b.getGeneSymbol())
      && a.getDepth() == b.getDepth();
  }

  private static boolean isSameMap(final Map<String, TreeMap<Integer, Interval>> a,
      final Map<String, TreeMap<Integer, Interval>> b) {
    if (a == null || b == null) {
      return false;
    }
    if (a.size() != b.size()) {
      return false;
    }
    for (final String chr : a.keySet()) {
      final TreeMap<Integer, Interval> atm = a.get(chr);
      final TreeMap<Integer, Interval> btm = b.get(chr);
      if (atm.size() != btm.size()) {
        return false;
      }
      for (final int start : atm.keySet()) {
        if (!isSameInterval(atm.get(start), btm.get(start))) {
          return false;
        }
      }
    }
    return true;
  }

  private static Stream<Arguments> toIntListArgs() {
    return Stream.of(
        Arguments.of(Arrays.asList(16400226), "16400226,"),
        Arguments.of(Arrays.asList(592235, 595220), "592235,595220,"));
  }

  @ParameterizedTest
  @MethodSource("toIntListArgs")
  void toIntListTest(final List<Integer> expected, final String input)
      throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
    final Method toIntList = GeneExons.class.getDeclaredMethod("toIntList", String.class);
    toIntList.setAccessible(true);
    assertIterableEquals(expected, (List<Integer>)toIntList.invoke(GeneExons.class, input));
  }

  @Test
  void loadmapTest() {
    final GeneExons ge = new GeneExons(refGene);
    assertEquals(expectedCounter, ge.getCounter());
    assertTrue(isSameMap(expectedGeneMap, ge.getGeneMap()));
    assertTrue(isSameMap(expectedMap, ge.getMap()));
  }
}
