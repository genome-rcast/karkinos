package jp.ac.utokyo.rcast.karkinos.wavelet;

import static java.lang.Double.NaN;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import jp.ac.utokyo.rcast.karkinos.cntwavelet.CTWaveletBean;
import jp.ac.utokyo.rcast.karkinos.cntwavelet.GaussianWavelet;
import jp.ac.utokyo.rcast.karkinos.exec.CNVInfo;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import manifold.ext.rt.api.Jailbreak;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class EMMethodTest {
    static Stream<Arguments> calcParameters() {
        double[][] deNoiseValues1 = new double[1][4000];
        IntStream.range(0, 4000)
                .forEach(i -> deNoiseValues1[0][i] = 0.001 * i);
        int[] expSignalCount1 = new int[4000];
        Arrays.fill(expSignalCount1, 1, 3999, 1);
        double[] expMa1 = new double[4000];
        IntStream.range(0, 6)
                .forEach(i -> expMa1[i] = (double)(5 + i) / (6 + i));
        Arrays.fill(expMa1, 6, 3994, 1.0);
        IntStream.range(0, 6)
                .forEach(i -> expMa1[3994+i] = (double)(10 - i) / (11 - i));

        double[][] deNoiseValues2 = new double[1][100];
        IntStream.range(0, 100)
                .forEach(i -> deNoiseValues2[0][i] = 0.001 * i);
        int[] expSignalCount2 = new int[4000];
        Arrays.fill(expSignalCount2, 1, 99, 1);
        double[] expMa2 = new double[4000];
        IntStream.range(0, 6)
                .forEach(i -> expMa2[i] = (double)(5 + i) / (6 + i));
        Arrays.fill(expMa2, 6, 94, 1.0);
        IntStream.range(0, 10)
                .forEach(i -> expMa2[94+i] = (double)(10 - i) / 11);

        return Stream.of(
                //deNoiseValues, baselineLOHEstimate,
                // expSignalCount, expMa, expPeaksXYs
                Arguments.of(
                        deNoiseValues1, 0.0, expSignalCount1, expMa1,
                        new double[][] {
                                {0.223, 0.04258915030662313},
                                {1.961, 0.04258915030662828},
                                {3.779, 0.04258915030662762},
                                {3.971, 0.03572781943467793}}),
                Arguments.of(
                        deNoiseValues2, 0.0, expSignalCount2, expMa2,
                        new double[][] {}));
    }

    @ParameterizedTest
    @MethodSource("calcParameters")
    void calcTest(final double[][] deNoiseValues,
                  final double baselineLOHEstimate,
                  final int[] expSignalCount,
                  final double[] expMa,
                  final double[][] expPeaksXYs) throws IOException {
        @Jailbreak final DataSet dataset = new DataSet(null);
        dataset.capintervalinit = true;
        dataset.capinterval = IntStream.range(0, deNoiseValues.length).mapToObj(i->
                IntStream.range(0, deNoiseValues[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, 0.0));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        int chr = 0;
        for (List<WaveletIF> l : dataset.capinterval) {
            int cnt = 0;
            for (WaveletIF wi : l) {
                wi.setDenioseValue(deNoiseValues[chr][cnt]);
                cnt++;
            }
            chr++;
        }
        final PeaksInfo pi = EMMethod.calc(dataset, baselineLOHEstimate);
        final double delta = 1e-6;

        CTWaveletBean bean = GaussianWavelet.getWaveletTransform(expMa);
        final double[] expPeakSignals = bean.getData();
        final double peakV = bean.getMostfittingvariance();

        assertArrayEquals(expSignalCount, pi.getSignalcount());
        assertArrayEquals(expMa, pi.getMa(), delta);
        assertArrayEquals(expPeakSignals, pi.getPeaksignals(), delta);

        final List<Peak> pl = pi.getPeaklist();
        assertEquals(expPeaksXYs.length, pl.size());

        double sumY = 0;
        for (double [] xy : expPeaksXYs) {
            sumY = sumY + xy[1];
        }
        int cnt = 0;
        for (Peak p : pl ) {
            assertEquals(expPeaksXYs[cnt][1], p.getY());
            assertEquals(expPeaksXYs[cnt][0], p.getU());
            assertEquals(peakV, p.getV());
            assertEquals(expPeaksXYs[cnt][1] / sumY, p.getR(), delta);
            assertEquals(cnt, p.getPeakidx());
            cnt++;
        }
    }

    static Stream<Arguments> isPeakParameters() {
        return Stream.of(
                //x, b4b, b4, v, next, nextNext, expected
                Arguments.of(0.0, 1.0, 2.0, 3.0, 2.0, 1.0, true),
                Arguments.of(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, false),
                Arguments.of(0.0, 1.0, 2.0, 3.0, 2.0, 2.0, false),
                Arguments.of(0.0, 1.0, 2.0, 2.0, 4.0, 5.0, false),
                Arguments.of(0.0, 1.0, 6.0, 10.0, 12.0, 14.0, false),
                Arguments.of(0.0, 2.0, 6.0, 10.0, 12.0, 13.0, false),
                Arguments.of(0.0, 1.0, 6.0, 10.0, 12.0, 12.0, false),
                Arguments.of(0.0, 1.0, 6.0, 10.0, 12.0, 13.0, true),
                Arguments.of(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, false));
    }

    @ParameterizedTest
    @MethodSource("isPeakParameters")
    void isPeakTest(final double x,
                    final double b4b,
                    final double b4,
                    final double v,
                    final double next,
                    final double nextNext,
                    final boolean expected) {
        @Jailbreak EMMethod emm = null;
        assertEquals(expected, emm.isPeak(x, b4b, b4, v, next, nextNext));
    }

    static Stream<Arguments> aveParameters() {
        return Stream.of(
                //bin, i, j, expected
                Arguments.of(
                        new int[] {1, 2, 3, 4, 5},
                        0, 4, 3.0),
                Arguments.of(
                        new int[] {1, 2, 3, 4, 5},
                        -1, 3, 2.5),
                Arguments.of(
                        new int[] {1, 2, 3, 4, 5},
                        3, 5, 4.5),
                Arguments.of(
                        new int[] {1, 2, 3, 4, 5},
                        2, 3, 3.5),
                Arguments.of(
                        new int[] {1, 2, 3, 4, 5},
                        3, -1, NaN));
    }

    @ParameterizedTest
    @MethodSource("aveParameters")
    void aveTest(final int[] bin,
                 final int i,
                 final int j,
                 final double expected) {
        @Jailbreak EMMethod emm = null;
        assertEquals(expected, emm.ave(bin, i, j));
    }

    static Stream<Arguments> getIdxParameters() {
        return Stream.of(
                //d, expected
                Arguments.of(-0.01, 0),
                Arguments.of(0.0, 0),
                Arguments.of(0.001, 1),
                Arguments.of(4.0, 3999),
                Arguments.of(4.001, 3999),
                Arguments.of(3.998, 3998),
                Arguments.of(0.0015, 2),
                Arguments.of(0.0014, 1));
    }

    @ParameterizedTest
    @MethodSource("getIdxParameters")
    void getIdxTest(final double d, final int expected) {
        @Jailbreak EMMethod emm = null;
        assertEquals(expected, emm.getIdx(d));
    }
}
