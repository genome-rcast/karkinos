package jp.ac.utokyo.rcast.karkinos.wavelet;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import manifold.ext.rt.api.Jailbreak;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.CNVInfo;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;

import static java.lang.Double.NaN;
import static org.junit.jupiter.api.Assertions.*;

public class WaveletDenoizeTest {
    static Stream<Arguments> calcParameters() {
        return Stream.of(
                // normalCounts, tnRatios, expLoh, expDeNoiseValues, expCNs
                Arguments.of(
                        new int[][] {{1, 2, 3}, {4, 5, 6}},
                        new double[][] {{0.1, 0.2, 0.3}, {1.1, 1.2, 1.3}},
                        0.51,
                        new double[][] {{0.21765567128029326, 0.21765567128029326,
                                0.21765567128029326}, {1.206723294623865,
                                1.206723294623865, 1.206723294623865}},
                        new double[][] {{0.5, 0.5, 0.5}, {1.0, 1.0, 1.0}}),
                Arguments.of(
                        new int[][]{{1, 2, 3, 4, 5, 6}},
                        new double[][] {{2.1, 2.2, 2.3, 2.1, 2.2, 2.3}},
                        0.90,
                        new double[][] {{2.2109080498231073, 2.2109080498231073,
                                2.2109080498231073, 2.2109080498231073,
                                2.2109080498231073, 2.2109080498231073}},
                        new double[][] {{7.0, 7.0, 7.0, 7.0, 7.0, 7.0}}));
    }

    @ParameterizedTest
    @MethodSource("calcParameters")
    void calcTest(final int[][] normalCounts,
                  final double [][] tnRatios,
                  final double expLoh,
                  final double[][] expDeNoiseValues,
                  final double[][] expCNs) throws IOException {
        @Jailbreak final DataSet dataset = new DataSet(null);
        dataset.capintervalinit = true;
        dataset.capinterval = IntStream.range(0, tnRatios.length).mapToObj(i->
                IntStream.range(0, tnRatios[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(normalCounts[i][j],
                            0.0, 0.0, 1, 0.0, 0.0, tnRatios[i][j]));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        final double baselineLOHEstimate = WaveletDenoize.calc(dataset);

        final double delta = 1e-7;
        assertEquals(expLoh, baselineLOHEstimate, delta);
        final List<List<WaveletIF>> capInterval = dataset.getCapInterval();
        int chr = 0;
        for (final List<WaveletIF> l : capInterval) {
            assertArrayEquals(expDeNoiseValues[chr],
                    l.stream().mapToDouble(WaveletIF::getDenioseValue).toArray());
            assertArrayEquals(expCNs[chr],
                    l.stream().mapToDouble(WaveletIF::getCN).toArray());
            chr++;
        }
    }

    static Stream<Arguments> getDenoiseLevelParameters() {
        double[][] TnRatio1 = new double[1][64];
        Arrays.fill(TnRatio1[0], 2.0);

        double[][] TnRatio2 = new double[1][4096];
        IntStream.range(0, 4096 / 2)
                .filter(i -> i % 2 != 0)
                .forEach(i -> TnRatio2[0][i] = 0.05);
        IntStream.range(4096 / 2, 4096)
                .filter(i -> i % 2 != 0)
                .forEach(i -> TnRatio2[0][i] = 0.04);

        double[][] TnRatio3 = new double[1][4096];
        Arrays.fill(TnRatio3[0], 1.0);

        return Stream.of(
                // tnRatios, expected
                Arguments.of(TnRatio1, 8),
                Arguments.of(TnRatio2, 2),
                Arguments.of(TnRatio3, 1));
    }

    @ParameterizedTest
    @MethodSource("getDenoiseLevelParameters")
    void getDenoiseLevelTest(final double [][] tnRatios,
                             final int expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, tnRatios.length).mapToObj(i->
                IntStream.range(0, tnRatios[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, tnRatios[i][j]));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoize wd = null;
        assertEquals(expected, wd.getDenoiseLevel(plist));
    }

    static Stream<Arguments> getSamplingDataParameters() {
        double[][] TnRatio = new double[][] {
                {0.5, 0.7},
                {1.5, 2.6},
                {2.4}};

        return Stream.of(
                // tnRatios, samplingSize, expected
                Arguments.of(
                        TnRatio,
                        5,
                        new double[] {0.5, 0.7, 1.5, 2.6, 2.4}),
                Arguments.of(
                        TnRatio,
                        4,
                        new double[] {0.5, 0.7, 1.5, 2.6}),
                Arguments.of(
                        TnRatio,
                        3,
                        new double[] {0.5, 0.7, 1.5}),
                Arguments.of(
                        TnRatio,
                        6,
                        new double[] {0.5, 0.7, 1.5, 2.6, 2.4, 0.0}));
    }

    @ParameterizedTest
    @MethodSource("getSamplingDataParameters")
    void getSamplingDataTest(final double [][] tnRatios,
                             final int samplingSize,
                             final double[] expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, tnRatios.length).mapToObj(i->
                IntStream.range(0, tnRatios[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, tnRatios[i][j]));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoize wd = null;
        double[] sampling = wd.getSamplingData(plist, samplingSize);

        final double delta = 1e-7;
        assertEquals(samplingSize, sampling.length);
        assertArrayEquals(expected, sampling, delta);
    }

    static Stream<Arguments> downsamplingParameters() {
        return Stream.of(
                //testSampling, expected
                Arguments.of(
                        new double[] {0.1, 0.2, 0.3, 0.4,
                                0.5, 0.6, 0.7, 0.8, 0.9},
                        new double[] {0.15, 0.35, 0.55, 0.75}));
    }

    @ParameterizedTest
    @MethodSource("downsamplingParameters")
    void downsamplingTest(final double[] testSampling, final double[] expected) {
        @Jailbreak WaveletDenoize wd = null;
        final double[] downSampling = wd.downsampling(testSampling);

        final double delta = 1e-7;
        assertEquals(expected.length, downSampling.length);
        assertArrayEquals(expected, downSampling, delta);
    }

    static Stream<Arguments> getSDParameters() {
        return Stream.of(
                //double[] testSampling, double expected
                Arguments.of(new double[] {0.1, 0.1, 0.1}, 0.0),
                Arguments.of(new double[] {0.1, 0.05, 0.15}, 0.05));
    }

    @ParameterizedTest
    @MethodSource("getSDParameters")
    void getSDTest(final double[] testSampling, final double expected) {
        @Jailbreak WaveletDenoize wd = null;
        final SummaryStatistics ss = wd.getSD(testSampling);
        final double delta = 1e-7;
        assertEquals(expected, ss.getStandardDeviation(), delta);
    }

    static Stream<Arguments> setcnParameters() {
        return Stream.of(
                //tnRatio, deNoiseValue, baselineLOHEstimate, expected
                Arguments.of(
                        new double[] {0.5, 0.7, 1.5, 2.6, 2.4},
                        new double[] {0.9, 0.7, 1.5, 3.2, 2.4},
                        0.6,
                        new double[] {0.5, 0.5, 1.5, 2.5, 2.5}));
    }

    @ParameterizedTest
    @MethodSource("setcnParameters")
    void setcnTest(final double [] tnRatio,
                   final double[] deNoiseValue,
                   final double baselineLOHEstimate,
                   final double[] expected) {
        List<WaveletIF> list = IntStream.range(0, tnRatio.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, tnRatio[j]));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        setDeNoiseValue(list, deNoiseValue);

        @Jailbreak WaveletDenoize wd = null;
        wd.setcn(list, baselineLOHEstimate);

        assertArrayEquals(expected,
                list.stream().mapToDouble(WaveletIF::getCN).toArray());
    }

    static Stream<Arguments> lohEstimateParameters() {
        return Stream.of(
                //deNoiseValues, expected
                Arguments.of(
                        new double[][] {{1.577, 1.578, 0.45, 0.4, 0.48}},
                        0.44),
                Arguments.of(
                        new double[][] {{2.077, 2.078, 2.079, 0.48, 0.2}},
                        0.58),
                Arguments.of(
                        new double[][] {{1.577, 1.578, 0.65, 0.6, 0.68}},
                        0.53),
                Arguments.of(
                        new double[][] {{1.577, 1.578, 0.65, 0.64, 0.66}},
                        0.69),
                Arguments.of(
                        new double[][] {{0.39, 0.38, 0.37, 0.36, 0.35}},
                        0.90));
    }

    @ParameterizedTest
    @MethodSource("lohEstimateParameters")
    void lohEstimateTest(final double[][] deNoiseValues, final double expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, deNoiseValues.length).mapToObj(i->
                IntStream.range(0, deNoiseValues[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, 0.0));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        int chr = 0;
        for (List<WaveletIF> l : plist) {
            setDeNoiseValue(l, deNoiseValues[chr]);
            chr++;
        }

        @Jailbreak WaveletDenoize wd = null;
        final double delta = 1e-6;
        assertEquals(expected, wd.lohEstimate(plist), delta);
    }

    static Stream<Arguments> lohEstimate2Parameters() {
        return Stream.of(
                //deNoiseValues, startVal, endVal, include3n, expected
                Arguments.of(
                        new double[][] {{0.751, 0.751, 0.751, 0.751, 0.751}},
                        0.45f, 0.75f, false, 0.0),
                Arguments.of(
                        new double[][] {{0.5, 0.6, 0.5, 0.6, 0.5}},
                        0.45f, 0.75f, false, 0.55),
                Arguments.of(
                        new double[][] {{2.077, 2.078, 2.079, 2.076, 2.080}},
                        0.45f, 0.75f, false, 0.53),
                Arguments.of(
                        new double[][] {{2.077, 2.078, 2.079, 2.076, 2.081}},
                        0.45f, 0.75f, false, 0.64),
                Arguments.of(
                        new double[][] {{0.9, 0.83, 0.8, 0.75, 0.7}},
                        0.65f, 0.95f, false, 0.78));
    }

    @ParameterizedTest
    @MethodSource("lohEstimate2Parameters")
    void lohEstimateTest(final double[][] deNoiseValues,
                         final float startVal,
                         final float endVal,
                         final boolean include3n,
                         final double expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, deNoiseValues.length).mapToObj(i->
                IntStream.range(0, deNoiseValues[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, 0.0));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        int chr = 0;
        for (List<WaveletIF> l : plist) {
            setDeNoiseValue(l, deNoiseValues[chr]);
            chr++;
        }

        @Jailbreak WaveletDenoize wd = null;
        final double delta = 1e-6;
        assertEquals(
                expected,
                wd.lohEstimate(plist, startVal, endVal, include3n),
                delta);
    }

    static Stream<Arguments> getLowerPeakParameters() {
        double[][] deNoiseValue1 = new double[1][100];
        Arrays.fill(deNoiseValue1[0], 0, 96, 0.7);
        Arrays.fill(deNoiseValue1[0], 96, 100, 0.24);

        double[][] deNoiseValue2 = new double[1][100];
        Arrays.fill(deNoiseValue2[0], 0, 97, 0.7);
        Arrays.fill(deNoiseValue2[0], 97, 100, 0.24);

        List<double[]> l1 = Arrays.asList(
                new double[]{0.1, 0.9},
                new double[]{0.2, 0.8},
                new double[]{0.3, 0.7},
                new double[]{0.4, 0.6},
                new double[]{0.4, 0.5});
        List<double[]> l2 = Arrays.asList(
                new double[]{0.1, 0.9},
                new double[]{0.3, 0.5},
                new double[]{0.5, 0.7},
                new double[]{0.2, 0.1},
                new double[]{0.6, 0.2});
        List<double[]> l3 = Arrays.asList(
                new double[]{0.1, 0.9},
                new double[]{0.25, 0.5},
                new double[]{0.5, 0.7},
                new double[]{0.3, 0.1},
                new double[]{0.7, 0.2});
        List<double[]> l4 = Arrays.asList(
                new double[]{0.1, 0.9},
                new double[]{0.2, 0.5},
                new double[]{0.5, 0.7},
                new double[]{0.3, 0.1},
                new double[]{0.7, 0.2});
        List<double[]> l5 = Arrays.asList(
                new double[]{0.1, 0.9},
                new double[]{0.2, 0.5},
                new double[]{0.5, 0.7},
                new double[]{0.3, 0.1},
                new double[]{0.7, 0.2});

        return Stream.of(
                //l, l_low, deNoiseValues, expected
                Arguments.of(
                        l1, l1,
                        new double[][] {{0.0}},
                        -1.0f),
                Arguments.of(
                        l2, l2,
                        new double[][] {{0.0}},
                        0.2f),
                Arguments.of(
                        l3, l3,
                        new double[][] {{0.0}},
                        0.3f),
                Arguments.of(
                        l4, l4,
                        deNoiseValue1,
                        0.2f),
                Arguments.of(
                        l5, l5,
                        deNoiseValue2,
                        0.3f));
    }

    @ParameterizedTest
    @MethodSource("getLowerPeakParameters")
    void getLowerPeakTest(final List<double[]> l,
                          final List<double[]> l_low,
                          final double[][] deNoiseValues,
                          final float expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, deNoiseValues.length).mapToObj(i->
                IntStream.range(0, deNoiseValues[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, 0.0));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        int chr = 0;
        for (List<WaveletIF> lwi : plist) {
            setDeNoiseValue(lwi, deNoiseValues[chr]);
            chr++;
        }

        @Jailbreak WaveletDenoize wd = null;
        assertEquals(expected, wd.getLowerPeak(l, l_low, plist));
    }

    static Stream<Arguments> moreThanParcentLowAreaParameters() {
        double[][] deNoiseValues = new double[1][50];
        Arrays.fill(deNoiseValues[0], 0.7);
        deNoiseValues[0][0] = 0.4;

        return Stream.of(
                //f, firstPeak, deNoiseValues, expected
                Arguments.of(0.99, 0.71, deNoiseValues, true),
                Arguments.of(1.0, 0.71, deNoiseValues, false),
                Arguments.of(1.01, 0.71, deNoiseValues, false),
                Arguments.of(0.019, 0.70, deNoiseValues, true),
                Arguments.of(0.02, 0.70, deNoiseValues, false),
                Arguments.of(0.021, 0.70, deNoiseValues, false));
    }

    @ParameterizedTest
    @MethodSource("moreThanParcentLowAreaParameters")
    void moreThanParcentLowAreaTest(final double f,
                                    final double firstPeak,
                                    final double[][] deNoiseValues,
                                    final boolean expected) {
        List<List<WaveletIF>> plist = IntStream.range(0, deNoiseValues.length).mapToObj(i->
                IntStream.range(0, deNoiseValues[i].length).mapToObj(j->{
                    final CapInterval ci = new CapInterval("chr"+i, j, j+1, true);
                    ci.setCNVInfo(new CNVInfo(1,
                            0.0, 0.0, 1, 0.0, 0.0, 0.0));
                    return (WaveletIF)ci;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        int chr = 0;
        for (List<WaveletIF> l : plist) {
            setDeNoiseValue(l, deNoiseValues[chr]);
            chr++;
        }

        @Jailbreak WaveletDenoize wd = null;
        assertEquals(expected, wd.moreThanParcentLowArea(f, firstPeak, plist));
    }

    static Stream<Arguments> sortParameters() {
        double[] data1 = new double[] {0.1, 0.9};
        double[] data2 = new double[] {0.2, 0.5};
        double[] data3 = new double[] {0.3, 0.7};

        return Stream.of(
                //peaks, expected
                Arguments.of(
                        Arrays.asList(data1, data2, data3),
                        Arrays.asList(data2, data3, data1)));
    }

    @ParameterizedTest
    @MethodSource("sortParameters")
    void sortTest(final List<double[]> peaks, final List<double[]> expected) {
        @Jailbreak WaveletDenoize wd = null;
        wd.sort(peaks);
        assertEquals(expected, peaks);
    }

    static Stream<Arguments> getNearistLineParameters() {
        return Stream.of(
                //interval, val, include3n, expected
                Arguments.of(0.15, 1.22, true, 1.15),
                Arguments.of(0.15, 1.225, true, 1.3),
                Arguments.of(0.15, 1.23, true, 1.3),
                Arguments.of(0.15, 1.3, false, 1.3)
        );
    }

    @ParameterizedTest
    @MethodSource("getNearistLineParameters")
    void getNearistLineTest(final double interval,
                            final double val,
                            final boolean include3n,
                            final double expected) {
        @Jailbreak WaveletDenoize wd = null;
        assertEquals(expected, wd.getNearistLine(interval, val, include3n));
    }

    static Stream<Arguments> _calcMovingAverageParameters() {
        return Stream.of(
                //normalCount, tnRatio, deNoiseLevel, expected
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        1,
                        new double[] {
                                0.15857864376269054, 0.15857864376269054,
                                0.25505102572168215, 0.3535898384862246,
                                0.4527864045000421}),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        8,
                        new double[] {
                                0.3364805672917484, 0.3364805672917484,
                                0.3364805672917484, 0.3364805672917484,
                                0.3364805672917484}),
                Arguments.of(
                        new long[] {0},
                        new double[] {0.0},
                        8,
                        new double[] {0.0}));
    }

    @ParameterizedTest
    @MethodSource("_calcMovingAverageParameters")
    void _calcMovingAverageTest(final long[] normalCount,
                                final double[] tnRatio,
                                final int deNoiseLevel,
                                final double[] expected) {
        List<WaveletIF> list = IntStream.range(0, tnRatio.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(normalCount[j],
                    0.0, 0.0, 1, 0.0, 0.0, tnRatio[j]));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoize wd = null;
        wd._calcMovingAverage(list, deNoiseLevel);
        final double delta = 1e-7;
        assertArrayEquals(expected,
                list.stream().mapToDouble(WaveletIF::getDenioseValue).toArray(),
                delta);
    }

    static Stream<Arguments> getMAParameters() {
        return Stream.of(
                //normalCount, tnRatio, idx, size, expected
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        0, 0, NaN),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        0, 8, 0.3364805672917484),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        0, 2, 0.15857864376269054),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        1, 2, 0.15857864376269054),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        2, 2, 0.25505102572168215),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        3, 2, 0.3535898384862246),
                Arguments.of(
                        new long[] {1, 2, 3, 4, 5},
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        4, 2, 0.4527864045000421));
    }

    @ParameterizedTest
    @MethodSource("getMAParameters")
    void getMATest(final long[] normalCount,
                   final double[] tnRatio,
                   final int idx,
                   final int size,
                   final double expected) {
        List<WaveletIF> list = IntStream.range(0, tnRatio.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(normalCount[j],
                    0.0, 0.0, 1, 0.0, 0.0, tnRatio[j]));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoize wd = null;
        final double delta = 1e-7;
        assertEquals(expected, wd.getMA(idx, size, list), delta);
    }

    static Stream<Arguments> adjustBoundaryParameters() {
        return Stream.of(
                //tnRatio, cn, expected
                Arguments.of(
                        new double[] {0.1, 0.2, 0.3, 0.4, 0.5},
                        new double[] {0.16, 0.14, 0.3, 0.46, 0.45},
                        new double[] {0.14, 0.14, 0.3, 0.45, 0.45}));
    }

    @ParameterizedTest
    @MethodSource("adjustBoundaryParameters")
    void adjustBoundaryTest(final double[] tnRatio,
                            final double[] cn,
                            final double[] expected) {
        List<WaveletIF> list = IntStream.range(0, tnRatio.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, tnRatio[j]));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        setCN(list, cn);

        @Jailbreak WaveletDenoize wd = null;
        wd.adjustBoundary(list);
        assertArrayEquals(expected,
                list.stream().mapToDouble(WaveletIF::getCN).toArray());
    }

    static Stream<Arguments> boundCheckParameters() {
        return Stream.of(
                //tnRatio, cn, n, expected
                Arguments.of(
                        new double[] {0.5, 0.5, 0.5, 0.5, 0.5},
                        new double[] {0.6, 0.5, 0.6, 0.51, 0.6},
                        3,
                        new double[] {0.6, 0.5, 0.51, 0.51, 0.6}),
                Arguments.of(
                        new double[] {0.5, 0.5, 0.5, 0.5, 0.5},
                        new double[] {0.51, 0.52, 0.6, 0.5, 0.6},
                        1,
                        new double[] {0.51, 0.51, 0.51, 0.5, 0.6}));
    }

    @ParameterizedTest
    @MethodSource("boundCheckParameters")
    void boundCheckTest(final double[] tnRatio,
                        final double[] cn,
                        final int n,
                        final double[] expected) {
        List<WaveletIF> list = IntStream.range(0, tnRatio.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, tnRatio[j]));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        setCN(list, cn);

        @Jailbreak WaveletDenoize wd = null;
        wd.boundCheck(list, n);
        assertArrayEquals(expected,
                list.stream().mapToDouble(WaveletIF::getCN).toArray());
    }

    static Stream<Arguments> stepParameters() {
        return Stream.of(
                //d, loh, expected
                Arguments.of(0.9, 0.6, 1.0),
                Arguments.of(0.8, 0.8, 0.5),
                Arguments.of(3.5, 0.9, 13.5));
    }

    @ParameterizedTest
    @MethodSource("stepParameters")
    void stepTest(final double d, final double loh, final double expected) {
        @Jailbreak WaveletDenoize wd = null;
        assertEquals(expected, wd.step(d, loh));
    }

    private static void setCN(List<WaveletIF> list, final double[] cn) {
        int cnt = 0;
        for (WaveletIF wi : list) {
            wi.setCN(cn[cnt]);
            cnt++;
        }
    }

    private static void setDeNoiseValue(List<WaveletIF> list,
                                        final double[] deNoiseValue) {
        int cnt = 0;
        for (WaveletIF wi : list) {
            wi.setDenioseValue(deNoiseValue[cnt]);
            cnt++;
        }
    }
}
