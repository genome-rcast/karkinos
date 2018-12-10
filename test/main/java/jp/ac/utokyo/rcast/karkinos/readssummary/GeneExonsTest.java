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
    put("chr15", new TreeMap<Integer, Interval>() {{
      final Interval iv1 = new Interval("chr15", 49448194, 49448213, "NM_001001556", "120");
      iv1.exonidx = 1;
      put(49448194, iv1);

      final Interval iv2 = new Interval("chr15", 49493359, 49493447, "NM_001001556", "120");
      iv2.exonidx = 2;
      put(49493359, iv2);

      final Interval iv3 = new Interval("chr15", 49509387, 49509510, "NM_001001556", "120");
      iv3.exonidx = 3;
      put(49509387, iv3);

      final Interval iv4 = new Interval("chr15", 49528048, 49528138, "NM_001001556", "120");
      iv4.exonidx = 4;
      put(49528048, iv4);

      final Interval iv5 = new Interval("chr15", 49531418, 49531564, "NM_001001556", "120");
      iv5.exonidx = 5;
      put(49531418, iv5);

      final Interval iv6 = new Interval("chr15", 49574184, 49574282, "NM_001001556", "120");
      iv6.exonidx = 6;
      put(49574184, iv6);

      final Interval iv7 = new Interval("chr15", 49575763, 49575915, "NM_001001556", "120");
      iv7.exonidx = 7;
      put(49575763, iv7);

      final Interval iv8 = new Interval("chr15", 49584524, 49584734, "NM_001001556", "120");
      iv8.exonidx = 8;
      put(49584524, iv8);

      final Interval iv9 = new Interval("chr15", 49611801, 49612002, "NM_001001556", "120");
      iv9.exonidx = 9;
      put(49611801, iv9);

      final Interval iv10 = new Interval("chr15", 49620149, 49620356, "NM_001001556", "120");
      iv10.exonidx = 10;
      put(49620149, iv10);
    }});
  }};

  private static boolean isSameInterval(final Interval a, final Interval b) {
    return a.getChr().equals(b.getChr())
      && a.getStart() == b.getStart()
      && a.getEnd() == b.getEnd()
      && a.getRefseqid().equals(b.getRefseqid())
      && a.getExonidx() == b.getExonidx()
      && a.getGeneSymbol().equals(b.getGeneSymbol());
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
