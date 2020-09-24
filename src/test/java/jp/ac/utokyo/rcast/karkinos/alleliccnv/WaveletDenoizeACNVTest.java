package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import manifold.ext.rt.api.Jailbreak;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static java.lang.Float.NaN;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class WaveletDenoizeACNVTest {
    static Stream<Arguments> calcParameters() {
        return Stream.of(
                // rows, gcAdjustedHs, gcAdjustedLs, expWtValHLs
                Arguments.of(
                        new float[][] {{1f, 1.5f}},
                        new float[][] {{1f, 3f}},
                        new float[][] {{1f, 2f}},
                        new float[][][] {{{NaN, 1.0f}, {NaN, 1.0f}}}),
                Arguments.of(
                        new float[][] {{0.9f, 0.9f, 0.9f, 0.9f, 0.9f, 0.9f}},
                        new float[][] {{0.9f, 0.9f, 0.9f, 0.9f, 0.9f, 0.9f}},
                        new float[][] {{0.9f, 0.9f, 0.9f, 0.9f, 0.9f, 0.9f}},
                        new float[][][] {{{1f, 1f}, {1f, 1f}, {1f, 1f},
                                {1f, 1f}, {1f, 1f}, {1f, 1f}}}),
                Arguments.of(
                        new float[][] {{0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f}},
                        new float[][] {{0.4f, 0.6f, 0.5f, 0.7f, 0.4f, 0.6f,
                                0.8f, 0.4f, 0.9f, 0.4f, 0.5f, 0.6f,
                                0.3f, 0.7f, 0.4f, 0.6f, 0.5f, 0.6f,
                                0.8f, 0.4f, 0.6f, 0.4f, 0.5f, 0.6f}},
                        new float[][] {{0.9f, 0.5f, 0.6f, 0.9f, 0.8f, 0.6f,
                                0.5f, 0.7f, 0.4f, 0.9f, 0.8f, 0.6f,
                                0.5f, 0.8f, 0.6f, 0.9f, 0.8f, 0.6f,
                                0.7f, 0.5f, 0.3f, 0.9f, 0.8f, 0.6f}},
                        new float[][][] {{{1f, 0.9f}, {1f, 0.9f}, {1f, 0.9f},
                                {1f, 0.9f}, {1f, 0.9f}, {1.0582782f, 0.83327824f},
                                {1.0332782f, 0.88327825f}, {1.0332782f, 0.9082782f},
                                {1.0207782f, 0.87077826f}, {1.0082783f, 0.83327824f},
                                {1.0207782f, 0.8582782f}, {1f, 0.90000004f},
                                {1f, 0.9f}, {1f, 1f}, {1f, 0.93749994f},
                                {1.0082783f, 0.89577824f}, {1f, 0.9f},
                                {1.0207782f, 0.8582782f}, {1f, 0.8875f},
                                {1f, 0.9f}, {1f, 0.86249995f},
                                {1f, 0.86249995f}, {1f, 0.86249995f},
                                {1f, 0.86249995f}}}),
                Arguments.of(
                        new float[][] {{0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f}},
                        new float[][] {{0.1f, 0.1f, 1.6f, 1.6f, 1.6f, 1.6f,
                                1.6f, 1.6f, 1.6f, 1.6f, 1.6f, 1.6f,
                                0.01f, 1.6f, 1.6f, 1.6f, 1.6f, 1.6f,
                                1.6f, 1.6f, 1.6f, 1.6f, 1.6f, 1.6f}}, 
                        new float[][] {{0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f,
                                0.4f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f,
                                1.6f, 0.4f, 0.4f, 0.4f, 0.4f, 0.4f,
                                0.4f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f}},
                        new float[][][] {{{0.98866963f, 0.878152f},
                                {0.98866963f, 0.878152f}, {0.98866963f, 0.878152f},
                                {0.98866963f, 0.878152f}, {0.98866963f, 0.878152f},
                                {1.0761696f, 0.790652f}, {1.2636696f, 0.803152f},
                                {1.2636696f, 0.815652f}, {1.2636696f, 0.828152f},
                                {1.0649196f, 0.990652f}, {1.0617676f, 1.0f},
                                {1.0617676f, 1.0f}, {1.0617676f, 1.0f},
                                {1.0617676f, 1.0f}, {1.0617676f, 1.0f},
                                {1.0617676f, 1.0f}, {1.0649196f, 0.990652f},
                                {1.2636696f, 0.828152f}, {1.2636696f, 0.815652f},
                                {1.2636696f, 0.803152f}, {1.2636696f, 0.790652f},
                                {1.2636696f, 0.790652f}, {1.2636696f, 0.790652f},
                                {1.2636696f, 0.790652f}}}),
                Arguments.of(
                        new float[][] {{0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f}},
                        new float[][] {{1.4f, 1.5f, 1.3f, 1.1f, 1.4f, 1.6f,
                                1.2f, 1.4f, 0.9f, 1.4f, 1.5f, 1.6f,
                                1.3f, 0.7f, 1.4f, 1.6f, 1.5f, 1.6f,
                                0.8f, 1.4f, 1.6f, 1.4f, 1.5f, 1.6f}},
                        new float[][] {{0.9f, 0.5f, 0.6f, 0.9f, 0.8f, 0.6f,
                                0.5f, 0.7f, 0.4f, 0.9f, 0.8f, 0.6f,
                                0.5f, 0.8f, 0.6f, 0.9f, 0.8f, 0.6f,
                                0.7f, 0.5f, 0.3f, 0.9f, 0.8f, 0.6f}},
                        new float[][][] {{{1.0088658f, 0.9338659f},
                                {1.0088658f, 0.9338659f}, {1.0088658f, 0.9338659f},
                                {1.0088658f, 0.9338659f}, {1.0088658f, 0.9338659f},
                                {1.0f, 0.9250001f}, {1.0f, 0.9875f}, {1.0f, 0.9875f},
                                {1.0213659f, 0.9088659f}, {1.0088658f, 0.8713659f},
                                {0.9963659f, 0.9963659f}, {1f, 0.9875001f},
                                {1f, 0.9875001f}, {1.0213659f, 0.9838659f},
                                {1.0463659f, 0.9463659f}, {1.0f, 0.975f}, {1f, 0.9875f},
                                {1.0f, 0.92499995f}, {1.0588659f, 0.9088659f},
                                {1.0713658f, 0.9338659f}, {1.0713658f, 0.8963659f},
                                {1.0713658f, 0.8963659f}, {1.0713658f, 0.8963659f},
                                {1.0713658f, 0.8963659f}}}),
                Arguments.of(
                        new float[][] {{0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,
                                0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f}},
                        new float[][] {{1.4f, 1.5f, 1.3f, 1.1f, 1.4f, 1.6f,
                                1.2f, 1.4f, 0.9f, 1.4f, 1.5f, 1.6f,
                                1.3f, 0.7f, 1.4f, 1.6f, 1.5f, 1.6f,
                                0.8f, 1.4f, 1.6f, 1.4f, 1.5f, 1.6f}},
                        new float[][] {{1.4f, 1.5f, 1.3f, 1.1f, 1.4f, 1.6f,
                                1.2f, 1.4f, 0.9f, 1.4f, 1.5f, 1.6f,
                                1.3f, 0.7f, 1.4f, 1.6f, 1.5f, 1.6f,
                                0.8f, 1.4f, 1.6f, 1.4f, 1.5f, 1.6f}},
                        new float[][][] {{{1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f},
                                {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f},
                                {1f, 1f}, {0.9963659f, 0.9963659f}, {1f, 1f},
                                {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f},
                                {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f}, {1f, 1f},
                                {1f, 1f}, {1f, 1f}}}));
    }

    @ParameterizedTest
    @MethodSource("calcParameters")
    void calcTest(final float[][] rows,
                  final float[][] gcAdjustedHs,
                  final float[][] gcAdjustedLs,
                  final float[][][] expWtValHLs) throws IOException {
        List<List<SNVHolderPlusACnv>> pList = IntStream.range(0, rows.length).mapToObj(i ->
                IntStream.range(0, rows[i].length).mapToObj(j -> {
                    ACNVInfoBean higherA = new ACNVInfoBean();
                    ACNVInfoBean lowerA = new ACNVInfoBean();
                    higherA.gcadjusted = gcAdjustedHs[i][j];
                    lowerA.gcadjusted = gcAdjustedLs[i][j];
                    lowerA.row = rows[i][j];
                    PileUPResult pur = new PileUPResult();
                    SNVHolder snv = new SNVHolder();
                    snv.setNormal(pur);
                    snv.setTumor(pur);
                    SNVHolderPlusACnv d = new SNVHolderPlusACnv(snv, 0d);
                    d.highera = higherA;
                    d.lowera = lowerA;
                    return d;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());

        WaveletDenoizeACNV.calc(pList);
        final double delta = 1e-6;

        int i = 0;
        for (final List<SNVHolderPlusACnv> l : pList) {
            for (int j = 0; j < expWtValHLs[i].length; j++) {
                assertEquals(expWtValHLs[i][j][0], l.get(j).highera.wtval, delta);
                assertEquals(expWtValHLs[i][j][1], l.get(j).lowera.wtval, delta);
            }
            i++;
        }
    }

    static Stream<Arguments> getDenoiseLevelParameters() {
        float[][] rows1 = new float[1][64];
        Arrays.fill(rows1[0], 2.0f);

        float[][] rows2 = new float[1][4096];
        IntStream.range(0, 4096 / 2)
                .filter(i -> i % 2 != 0)
                .forEach(i -> rows2[0][i] = 0.05f);
        IntStream.range(4096 / 2, 4096)
                .filter(i -> i % 2 != 0)
                .forEach(i -> rows2[0][i] = 0.04f);

        float[][] rows3 = new float[1][4096];
        Arrays.fill(rows3[0], 1.0f);

        return Stream.of(
                // tnRatios, expected
                Arguments.of(rows1, 8),
                Arguments.of(rows2, 2),
                Arguments.of(rows3, 1));
    }

    @ParameterizedTest
    @MethodSource("getDenoiseLevelParameters")
    void getDenoiseLevelTest(final float[][] rows,
                             final int expected) {
        List<List<SNVHolderPlusACnv>> pList = Arrays.stream(rows)
                .map(row -> IntStream.range(0, row.length).mapToObj(j -> {
                    ACNVInfoBean lowerA = new ACNVInfoBean();
                    lowerA.row = row[j];
                    PileUPResult pur = new PileUPResult();
                    SNVHolder snv = new SNVHolder();
                    snv.setNormal(pur);
                    snv.setTumor(pur);
                    SNVHolderPlusACnv d = new SNVHolderPlusACnv(snv, 0d);
                    d.lowera = lowerA;
                    return d;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        assertEquals(expected, wdACnv.getDenoiseLevel(pList));
    }

    static Stream<Arguments> getSamplingDataParameters() {
        float[][] rows = new float[][] {
                {0.5f, 0.7f},
                {1.5f, 2.6f},
                {2.4f}};

        return Stream.of(
                // tnRatios, samplingSize, expected
                Arguments.of(
                        rows,
                        5,
                        new double[] {0.5, 0.7, 1.5, 2.6, 2.4}),
                Arguments.of(
                        rows,
                        4,
                        new double[] {0.5, 0.7, 1.5, 2.6}),
                Arguments.of(
                        rows,
                        3,
                        new double[] {0.5, 0.7, 1.5}),
                Arguments.of(
                        rows,
                        6,
                        new double[] {0.5, 0.7, 1.5, 2.6, 2.4, 0.0}));
    }

    @ParameterizedTest
    @MethodSource("getSamplingDataParameters")
    void getSamplingDataTest(final float[][] rows,
                             final int samplingSize,
                             final double[] expected) {
        List<List<SNVHolderPlusACnv>> pList = Arrays.stream(rows)
                .map(row -> IntStream.range(0, row.length).mapToObj(j -> {
                    ACNVInfoBean lowerA = new ACNVInfoBean();
                    lowerA.row = row[j];
                    PileUPResult pur = new PileUPResult();
                    SNVHolder snv = new SNVHolder();
                    snv.setNormal(pur);
                    snv.setTumor(pur);
                    SNVHolderPlusACnv d = new SNVHolderPlusACnv(snv, 0d);
                    d.lowera = lowerA;
                    return d;})
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());

        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        final double[] sampling = wdACnv.getSamplingData(pList, samplingSize);

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
        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        final double[] downSampling = wdACnv.downsampling(testSampling);

        final double delta = 1e-7;
        assertEquals(expected.length, downSampling.length);
        assertArrayEquals(expected, downSampling, delta);
    }

    static Stream<Arguments> getSDParameters() {
        return Stream.of(
                //testSampling, expected
                Arguments.of(new double[] {0.1, 0.1, 0.1}, 0.0),
                Arguments.of(new double[] {0.1, 0.05, 0.15}, 0.05));
    }

    @ParameterizedTest
    @MethodSource("getSDParameters")
    void getSDTest(final double[] testSampling, final double expected) {
        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        final SummaryStatistics ss = wdACnv.getSD(testSampling);
        final double delta = 1e-7;
        assertEquals(expected, ss.getStandardDeviation(), delta);
    }

    static Stream<Arguments> _calcMovingAverageParameters() {
        return Stream.of(
                //gcAdjustedHL, deNoiseLevel, expWtValHL
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f},
                                {5f, 6f}, {7f, 8f}, {9f, 10f}},
                        1,
                        new float[][] {{2.0f, 3.0f}, {2.0f, 3.0f},
                                {4.0f, 5.0f}, {6.0f, 7.0f}, {8.0f, 9.0f}}),
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f},
                                {5f, 6f}, {7f, 8f}, {9f, 10f}},
                        4,
                        new float[][] {{5.0f, 6.0f}, {5.0f, 6.0f}, 
                                {5.0f, 6.0f}, {5.0f, 6.0f}, {5.0f, 6.0f}}));
    }

    @ParameterizedTest
    @MethodSource("_calcMovingAverageParameters")
    void _calcMovingAverageTest(final float[][] gcAdjustedHL,
                                final int deNoiseLevel,
                                final float[][] expWtValHL) throws IOException {
        List<SNVHolderPlusACnv> dList = Arrays.stream(gcAdjustedHL).map(hl -> {
            ACNVInfoBean higherA = new ACNVInfoBean();
            ACNVInfoBean lowerA = new ACNVInfoBean();
            higherA.gcadjusted = hl[0];
            lowerA.gcadjusted = hl[1];
            PileUPResult pur = new PileUPResult();
            SNVHolder snv = new SNVHolder();
            snv.setNormal(pur);
            snv.setTumor(pur);
            SNVHolderPlusACnv d = new SNVHolderPlusACnv(snv, 0d);
            d.highera = higherA;
            d.lowera = lowerA;
            return d;}).collect(Collectors.toList());

        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        wdACnv._calcMovingAverage(dList, deNoiseLevel);
        for (int i = 0; i < expWtValHL.length; i++) {
            assertEquals(expWtValHL[i][0], dList.get(i).highera.wtval);
            assertEquals(expWtValHL[i][1], dList.get(i).lowera.wtval);
        }
    }

    static Stream<Arguments> getMAParameters() {
        return Stream.of(
                //gcAdjustedHL, idx, size, expected
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f}},
                        0, 0,
                        new double[] {0d, 0d}),
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f}, {5f, 6f}, {7f, 8f}},
                        0, 2,
                        new double[] {2.0, 3.0}),
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f}, {5f, 6f}, {7f, 8f}},
                        3, 2,
                        new double[] {6.0, 7.0}),
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f},
                                {5f, 6f}, {7f, 8f}, {9f, 10f}},
                        2, 6,
                        new double[] {5.0, 6.0}),
                Arguments.of(
                        new float[][] {{1f, 2f}, {3f, 4f},
                                {5f, 6f}, {7f, 8f}, {9f, 10f}},
                        2, 3,
                        new double[] {4.0, 5.0}));
    }

    @ParameterizedTest
    @MethodSource("getMAParameters")
    void getMATest(final float[][] gcAdjustedHL,
                   final int idx,
                   final int size,
                   final double[] expected) {
        List<SNVHolderPlusACnv> dList = Arrays.stream(gcAdjustedHL).map(hl -> {
            ACNVInfoBean higherA = new ACNVInfoBean();
            ACNVInfoBean lowerA = new ACNVInfoBean();
            higherA.gcadjusted = hl[0];
            lowerA.gcadjusted = hl[1];
            PileUPResult pur = new PileUPResult();
            SNVHolder snv = new SNVHolder();
            snv.setNormal(pur);
            snv.setTumor(pur);
            SNVHolderPlusACnv d = new SNVHolderPlusACnv(snv, 0d);
            d.highera = higherA;
            d.lowera = lowerA;
            return d;}).collect(Collectors.toList());

        @Jailbreak WaveletDenoizeACNV wdACnv = null;
        final double delta = 1e-6;
        assertArrayEquals(expected, wdACnv.getMA(idx, size, dList), delta);
    }
}
