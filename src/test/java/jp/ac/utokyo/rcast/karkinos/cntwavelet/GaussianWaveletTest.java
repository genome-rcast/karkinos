package jp.ac.utokyo.rcast.karkinos.cntwavelet;

import manifold.ext.rt.api.Jailbreak;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class GaussianWaveletTest {
    static Stream<Arguments> getWaveletTransformParameters() {
        double[] data1 = new double[4000];
        IntStream.range(0, 4000)
                .forEach(i -> data1[i] = 0.001 * i);
        double[] data2 = new double[4000];
        IntStream.range(0, 4000)
                .forEach(i -> data2[3999 - i] = 0.001 * i);
        double[] data3 = new double[4000];
        IntStream.range(0, 4000)
                .forEach(i -> data3[i] = 1.5);

        return Stream.of(
                //data, expSd
                Arguments.of(data1, 0.029),
                Arguments.of(data2, 0.029),
                Arguments.of(data3, 0.029));
    }

    @ParameterizedTest
    @MethodSource("getWaveletTransformParameters")
    void getWaveletTransformTest(final double[] data, final double expSd) {
        double sum = Arrays.stream(data).sum();
        double[] expData = new double[4000];
        IntStream.range(0, 4000)
                .forEach(i -> expData[i] = data[i] / sum);

        @Jailbreak GaussianWavelet gw = null;
        CTWaveletBean bean = GaussianWavelet.getWaveletTransform(data);

        final double delta = 1e-6;
        assertEquals(expSd, bean.getSd(), delta);
        assertEquals(Math.pow(expSd, 2), bean.getMostfittingvariance(), delta);
        assertEquals(
                gw.tf(data, expSd, 0),
                bean.getData()[0],
                delta);
        assertEquals(
                gw.tf(data, expSd, 3998),
                bean.getData()[3998],
                delta);
        assertArrayEquals(expData, data, delta);
    }

    static Stream<Arguments> tf2Parameters() {
        double[] data1 = new double[4000];
        IntStream.range(0, 4000)
                .forEach(i -> data1[i] = 0.001 * i);

        return Stream.of(
                //data, sd, n, expected
                Arguments.of(
                        data1,
                        0.01, 305, 30.5),
                Arguments.of(
                        data1,
                        0.01, 405, 40.5),
                Arguments.of(
                        data1,
                        0.02, 305, 43.133513652379364),
                Arguments.of(
                        data1,
                        0.02, 405, 57.27564927611027),
                Arguments.of(
                        new double[] {1.0},
                        0.2, -300, 0.2896101481912661),
                Arguments.of(
                        new double[] {1.0},
                        0.2, 3701, 0.0));
    }

    @ParameterizedTest
    @MethodSource("tf2Parameters")
    void tfTest(final double[] data,
                final double sd,
                final int n,
                final double expected) {
        @Jailbreak GaussianWavelet gw = null;
        final double delta = 1e-6;
        assertEquals(expected, gw.tf(data, sd, n), delta);
    }

    static Stream<Arguments> getGaussianParameters() {
        return Stream.of(
                //sd, diff, expected
                Arguments.of(2.0, 0.0, 0.28209479177387814),
                Arguments.of(0.2, 0.3, 0.2896101481912661),
                Arguments.of(0.2, 0.4, 0.12072747129440328),
                Arguments.of(0.3, 0.4, 0.2994400585271616),
                Arguments.of(0.001, 0.001, 7.651786165616441),
                Arguments.of(0.001, 0.003, 0.14014735226324268),
                Arguments.of(0.001, 0.004, 0.004232083331915876),
                Arguments.of(0.002, 0.003, 2.8961014819126607));
    }

    @ParameterizedTest
    @MethodSource("getGaussianParameters")
    void getGaussianTest(final double sd,
                         final double diff,
                         final double expected) {
        @Jailbreak GaussianWavelet gw = null;
        final double delta = 1e-6;
        assertEquals(expected, gw.getGaussian(sd, diff), delta);
    }

    static Stream<Arguments> pow2Parameters() {
        return Stream.of(
                //d, expected
                Arguments.of(2.0, 4),
                Arguments.of(0.2, 0.04),
                Arguments.of(1.414, 1.999396));
    }

    @ParameterizedTest
    @MethodSource("pow2Parameters")
    void pow2Test(final double d, final double expected) {
        final double delta = 1e-6;
        assertEquals(expected, GaussianWavelet.pow2(d), delta);
    }

    static Stream<Arguments> sqrtParameters() {
        return Stream.of(
                //d, expected
                Arguments.of(4.0, 2),
                Arguments.of(2.0, 1.41421356),
                Arguments.of(3.0, 1.7320508));
    }

    @ParameterizedTest
    @MethodSource("sqrtParameters")
    void sqrtTest(final double d, final double expected) {
        final double delta = 1e-6;
        assertEquals(expected, GaussianWavelet.sqrt(d), delta);
    }
}
