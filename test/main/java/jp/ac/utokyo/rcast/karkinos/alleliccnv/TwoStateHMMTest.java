package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.provider.MethodSource;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.List;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CNVInfo;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

// TODO: This test class only aims at testing CopyNumberInterval and CapInterval.
// We should add more tests in the future.
public class TwoStateHMMTest {
  @Test
  void getCoreVarianceTest() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
    final Method getCoreVariance = TwoStateHMM.class.getDeclaredMethod(
        "getCoreVariance", List.class, CopyNumberInterval.class);
    getCoreVariance.setAccessible(true);

    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {{
      setStart(4);
      setEnd(6);
    }};

    //      1234567890
    // cni     ===
    // ci1  ---
    // ci2  ---
    {
      CapInterval ci1 = new CapInterval("chr1", 1, 3, false);
      CNVInfo cnv1 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv1.setDenioseValue(1.0);
      ci1.setCNVInfo(cnv1);
      CapInterval ci2 = new CapInterval("chr1", 1, 3, false);
      CNVInfo cnv2 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv2.setDenioseValue(10.0);
      ci2.setCNVInfo(cnv2);
      final List<WaveletIF> neighbors = Arrays.asList(ci1, ci2);
      final SummaryStatistics ss = (SummaryStatistics)getCoreVariance.invoke(
          TwoStateHMM.class, neighbors, cni);

      SummaryStatistics expected = new SummaryStatistics();
      expected.addValue(1.0);
      expected.addValue(10.0);
      assertTrue(ss.equals(expected));
    }

    //      1234567890
    // cni     ===
    // ci1  ---
    // ci2   ---
    {
      CapInterval ci1 = new CapInterval("chr1", 1, 3, false);
      CNVInfo cnv1 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv1.setDenioseValue(1.0);
      ci1.setCNVInfo(cnv1);
      CapInterval ci2 = new CapInterval("chr1", 2, 4, false);
      CNVInfo cnv2 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv2.setDenioseValue(10.0);
      ci2.setCNVInfo(cnv2);
      final List<WaveletIF> neighbors = Arrays.asList(ci1, ci2);
      final SummaryStatistics ss = (SummaryStatistics)getCoreVariance.invoke(
          TwoStateHMM.class, neighbors, cni);

      SummaryStatistics expected = new SummaryStatistics();
      expected.addValue(10.0);
      assertTrue(ss.equals(expected));
    }

    //      1234567890
    // cni     ===
    // ci1  ---
    // ci2       ---
    {
      CapInterval ci1 = new CapInterval("chr1", 1, 3, false);
      CNVInfo cnv1 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv1.setDenioseValue(1.0);
      ci1.setCNVInfo(cnv1);
      CapInterval ci2 = new CapInterval("chr1", 6, 8, false);
      CNVInfo cnv2 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv2.setDenioseValue(10.0);
      ci2.setCNVInfo(cnv2);
      final List<WaveletIF> neighbors = Arrays.asList(ci1, ci2);
      final SummaryStatistics ss = (SummaryStatistics)getCoreVariance.invoke(
          TwoStateHMM.class, neighbors, cni);

      SummaryStatistics expected = new SummaryStatistics();
      expected.addValue(10.0);
      assertTrue(ss.equals(expected));
    }

    //      1234567890
    // cni     ===
    // ci1  ---
    // ci2        ---
    {
      CapInterval ci1 = new CapInterval("chr1", 1, 3, false);
      CNVInfo cnv1 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv1.setDenioseValue(1.0);
      ci1.setCNVInfo(cnv1);
      CapInterval ci2 = new CapInterval("chr1", 7, 9, false);
      CNVInfo cnv2 = new CNVInfo(0L, 0.0, 0.0, 0L, 0.0, 0.0, 0.0);
      cnv2.setDenioseValue(10.0);
      ci2.setCNVInfo(cnv2);
      final List<WaveletIF> neighbors = Arrays.asList(ci1, ci2);
      final SummaryStatistics ss = (SummaryStatistics)getCoreVariance.invoke(
          TwoStateHMM.class, neighbors, cni);

      SummaryStatistics expected = new SummaryStatistics();
      expected.addValue(1.0);
      expected.addValue(10.0);
      assertTrue(ss.equals(expected));
    }
  }
}
