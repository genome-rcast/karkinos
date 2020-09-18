package jp.ac.utokyo.rcast.karkinos.hmm;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import jp.ac.utokyo.rcast.karkinos.exec.CNVInfo;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.wavelet.Peak;
import jp.ac.utokyo.rcast.karkinos.wavelet.PeaksInfo;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;
import manifold.ext.rt.api.Jailbreak;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class HMMCNVAnalysisFromEMTest {
    static Stream<Arguments> calcParameters() {
        return Stream.of(
                //deNoiseValues, peakUVs, expPeakIdx, expHmmValue
                Arguments.of(
                        new double[][] {{1.0, 2.0, 3.0, 4.0, 5.0}},
                        new double[][] {},
                        new int[][] {{0, 0, 0, 0, 0}},
                        new double[][] {{1.0, 1.0, 1.0, 1.0, 1.0}}),
                Arguments.of(
                        new double[][] {{1.0, 2.0, 3.0, 4.0, 5.0}},
                        new double[][] {
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0}},
                        new int[][] {{0, 0, 0, 0, 0}},
                        new double[][] {{1.0, 1.0, 1.0, 1.0, 1.0}}),
                Arguments.of(
                        new double[][] {{1.0, 2.0, 3.0, 4.0, 5.0}},
                        new double[][] {
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {3.0, 4.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {20.0, 20.0}},
                        new int[][] {{1, 1, 1, 1, 1}},
                        new double[][] {{0.0, 0.0, 0.0, 0.0, 0.0}}),
                Arguments.of(
                        new double[][] {{0.5, 0.3, 0.2, 0.8, 0.6}},
                        new double[][] {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}},
                        new int[][] {{0, 0, 0, 0, 0}},
                        new double[][] {{0.0, 0.0, 0.0, 0.0, 0.0}}),
                Arguments.of(
                        new double[][] {{5.5, 5.8, 5.4, 5.2, 5.6}},
                        new double[][] {
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {20.0, 20.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {20.0, 20.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {20.0, 20.0},
                                {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {20.0, 20.0}},
                        new int[][] {{2, 2, 2, 2, 2}},
                        new double[][] {{0.0, 0.0, 0.0, 0.0, 0.0}}));
    }

    @ParameterizedTest
    @MethodSource("calcParameters")
    void calcTest(final double[][] deNoiseValues,
                  final double[][] peakUVs,
                  final int[][] expPeakIdx,
                  final double[][] expHmmValue) {
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
        List<Peak> peakList = new ArrayList<>();
        int cnt = 0;
        for (double[] uv : peakUVs) {
            Peak p = new Peak(0.0, 0.0);
            p.setU(uv[0]);
            p.setV(uv[1]);
            p.setPeakidx(cnt);
            peakList.add(p);
            cnt++;
        }
        PeaksInfo pi = new PeaksInfo();
        pi.setPeaklist(peakList);

        HMMCNVAnalysisFromEM.calc(dataset, pi);

        final List<List<WaveletIF>> capInterval = dataset.getCapInterval();
        chr = 0;
        for (final List<WaveletIF> l : capInterval) {
            assertArrayEquals(expPeakIdx[chr],
                    l.stream().mapToInt(i->((CapInterval)i).getPeakIdx()).toArray());
            assertArrayEquals(expHmmValue[chr],
                    l.stream().mapToDouble(WaveletIF::getHMMValue).toArray());
            chr++;
        }
    }

    static Stream<Arguments> setdefParameters() {
        return Stream.of(
                //num
                Arguments.of(1),
                Arguments.of(3));
    }

    @ParameterizedTest
    @MethodSource("setdefParameters")
    void setdefTest(final int num) {
        List<WaveletIF> list = IntStream.range(0, num).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, 0.0));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        hmmCnv.setdef(list);

        double[] expHmmValue = new double[num];
        Arrays.fill(expHmmValue, 1.0);

        assertArrayEquals(expHmmValue,
                list.stream().mapToDouble(WaveletIF::getHMMValue).toArray());
    }

    static Stream<Arguments> toPeakListIdxParameters() {
        return Stream.of(
                //peakIdx, expected
                Arguments.of(
                        new int[] {5, 4, 3, 2, 1},
                        new HashMap<Integer, Integer>() {
                            {
                                put(0, 5);
                                put(1, 4);
                                put(2, 3);
                                put(3, 2);
                                put(4, 1);
                            }}),
                Arguments.of(
                        new int[] {6, 9, 2, 5, 8},
                        new HashMap<Integer, Integer>() {
                            {
                                put(0, 6);
                                put(1, 9);
                                put(2, 2);
                                put(3, 5);
                                put(4, 8);
                            }}));
    }

    @ParameterizedTest
    @MethodSource("toPeakListIdxParameters")
    void toPeakListIdxTest(final int[] peakIdx,
                           final Map<Integer, Integer> expected) {
        List<Peak> peakList = new ArrayList<>();
        for (int d : peakIdx) {
            Peak p = new Peak(0.0, 0.0);
            p.setPeakidx(d);
            peakList.add(p);
        }

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        assertEquals(expected, hmmCnv.toPeakListIdx(peakList));
    }

    static Stream<Arguments> checkExsistanceParameters() {
        return Stream.of(
                //peakUVs, deNoiseValue, expPeakUV
                Arguments.of(
                        new double[][] {{1.0, 2.0}, {3.0, 4.0}},
                        new double[] {},
                        new double[][] {{1.0, 2.0}, {3.0, 4.0}}),
                Arguments.of(
                        new double[][] {{6.0, 1.0}, {12.0, 1.5}},
                        new double[] {12.0},
                        new double[][] {{12.0, 1.5}}),
                Arguments.of(
                        new double[][] {{6.0, 1.0}, {12.0, 1.0}},
                        new double[] {4.0},
                        new double[][] {{6.0, 1.0}}),
                Arguments.of(
                        new double[][] {
                                {6.0, 1.0}, {12.0, 1.0}, {5.0, 1.2}},
                        new double[] {4.0},
                        new double[][] {{6.0, 1.0}, {5.0, 1.2}}),
                Arguments.of(
                        new double[][] {{6.0, 1.0}, {12.0, 1.0}},
                        new double[] {4.0, 11.0},
                        new double[][] {{6.0, 1.0}, {12.0, 1.0}}));
    }

    @ParameterizedTest
    @MethodSource("checkExsistanceParameters")
    void checkExsistanceTest(final double[][] peakUVs,
                             final double[] deNoiseValue,
                             final double[][] expPeakUV) {
        List<Peak> peakList = new ArrayList<>();
        for (double[] uv : peakUVs) {
            Peak p = new Peak(0.0, 0.0);
            p.setU(uv[0]);
            p.setV(uv[1]);
            peakList.add(p);
        }
        List<WaveletIF> list = IntStream.range(0, deNoiseValue.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, 0.0));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        int cnt = 0;
        for (WaveletIF wi : list) {
            wi.setDenioseValue(deNoiseValue[cnt]);
            cnt++;
        }

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        final List<Peak> retPeakList = hmmCnv.checkExsistance(peakList, list);

        assertEquals(expPeakUV.length, retPeakList.size());
        assertArrayEquals(Arrays.stream(expPeakUV).mapToDouble(p -> p[0]).toArray(),
                retPeakList.stream().mapToDouble(Peak::getU).toArray());
        assertArrayEquals(Arrays.stream(expPeakUV).mapToDouble(p -> p[1]).toArray(),
                retPeakList.stream().mapToDouble(Peak::getV).toArray());
    }

    static Stream<Arguments> containParameters() {
        double[] deNoiseValue1 = new double[1000];
        Arrays.fill(deNoiseValue1, 10, 11, 1.0);
        double[] deNoiseValue2 = new double[1000];
        Arrays.fill(deNoiseValue2, 10, 12, 1.0);

        return Stream.of(
                //peakV, peakU, deNoiseValue, expected
                Arguments.of(
                        6.0, 4.0, deNoiseValue1, false),
                Arguments.of(
                        6.0, 4.0, deNoiseValue2, true),
                Arguments.of(
                        0.0, 0.0,
                        new double[] {},
                        false),
                Arguments.of(
                        6.0, 4.0,
                        new double[] {2.0},
                        true),
                Arguments.of(
                        6.0, 0.0,
                        new double[] {0.0, 12.0},
                        false));
    }

    @ParameterizedTest
    @MethodSource("containParameters")
    void containTest(final double peakU,
                     final double peakV,
                     final double[] deNoiseValue,
                     final boolean expected) {
        Peak p = new Peak(0.0, 0.0);
        p.setU(peakU);
        p.setV(peakV);
        List<WaveletIF> list = IntStream.range(0, deNoiseValue.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, 0.0));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        int cnt = 0;
        for (WaveletIF wi : list) {
            wi.setDenioseValue(deNoiseValue[cnt]);
            cnt++;
        }

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        assertEquals(expected, hmmCnv.contain(p, list));
    }

    static Stream<Arguments> getListParameters() {
        return Stream.of(
                //deNoiseValue, peakU, expPi, expPortion
                Arguments.of(
                        new double[] {-10.1, -12.0, -9.3},
                        new double[] {-10.0},
                        new double[] {-10.1, -10.0, -10.0}),
                Arguments.of(
                        new double[] {10.1, 12.0, 9.3},
                        new double[] {10.0},
                        new double[] {10.1, 10.0, 10.0}),
                Arguments.of(
                        new double[] {4.5, 13.5, 4.49, 13.51},
                        new double[] {13.0, 5.0},
                        new double[] {4.5, 13.5, 5.0, 13.0}));
    }

    @ParameterizedTest
    @MethodSource("getListParameters")
    void getListTest(final double[] deNoiseValue,
                     final double[] peakU,
                     final double[] expected) {
        List<WaveletIF> list = IntStream.range(0, deNoiseValue.length).mapToObj(j->{
            final CapInterval ci = new CapInterval("chr1", j, j+1, true);
            ci.setCNVInfo(new CNVInfo(1,
                    0.0, 0.0, 1, 0.0, 0.0, 0.0));
            return (WaveletIF)ci;})
                .collect(Collectors.toList());
        int cnt = 0;
        for (WaveletIF wi : list) {
            wi.setDenioseValue(deNoiseValue[cnt]);
            cnt++;
        }
        List<Peak> peakList = new ArrayList<>();
        for (double u : peakU) {
            Peak p = new Peak(0.0, 0.0);
            p.setU(u);
            peakList.add(p);
        }

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        List<ObservationReal> oList = hmmCnv.getList(list, peakList);
        assertEquals(deNoiseValue.length, oList.size());
        assertArrayEquals(expected,
                oList.stream().mapToDouble(o -> o.value).toArray());
    }

    static Stream<Arguments> getHMMParameters() {
        return Stream.of(
                //peakUVs, expPi, expPortion
                Arguments.of(
                        new double[][] {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}},
                        1d / 3, 0.1 / (3 - 1)));
    }

    @ParameterizedTest
    @MethodSource("getHMMParameters")
    void getHMMTest(final double[][] peakUVs,
                    final double expPi,
                    final double expPortion) {
        List<Peak> peakList = new ArrayList<>();
        for (double[] uv : peakUVs) {
            Peak p = new Peak(0.0, 0.0);
            p.setU(uv[0]);
            p.setV(uv[1]);
            peakList.add(p);
        }

        @Jailbreak HMMCNVAnalysisFromEM hmmCnv = null;
        Hmm<ObservationReal> hmm = hmmCnv.getHMM(peakList);

        for (int i = 0; i < peakList.size(); i++) {
            assertEquals(expPi, hmm.getPi(i));
            for (int j = 0; j < peakList.size(); j++) {
                if (i == j) {
                    assertEquals(0.9, hmm.getAij(i, j));
                } else {
                    assertEquals(expPortion, hmm.getAij(i, j));
                }
            }
        }
    }
}
