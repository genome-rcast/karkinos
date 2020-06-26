package jp.ac.utokyo.rcast.karkinos.bean;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static jp.ac.utokyo.rcast.karkinos.bean.BaitSampling.*;

import org.junit.Test;
import java.util.Arrays;

public class BailSamplingBeanTest {
    @Test
    public void ForwardPosTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(StartForward, 0, 2, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setForward(0, 1, false);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }

    @Test
    public void ForwardNegTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(StartReverse, bait_sample_length - 10, bait_sample_length, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setForward(bait_sample_length - 10, bait_sample_length - 1, true);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }

    @Test
    public void ReversePosTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(EndForward, 10, 21, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setReverse(10, 20, false);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }

    @Test
    public void ReverseNegTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(EndReverse, 35, 52, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setReverse(35, 51, true);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }

    @Test
    public void ForwardPosOverRangeTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(StartForward, 0, bait_sample_length, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setForward(-10, bait_sample_length + 10, false);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }

    @Test
    public void ReverseNegOverRangeTest() {
        int[] StartForward = new int[bait_sample_length];
        int[] StartReverse = new int[bait_sample_length];
        int[] EndForward = new int[bait_sample_length];
        int[] EndReverse = new int[bait_sample_length];

        Arrays.fill(EndReverse, 0, bait_sample_length, 1);
        BailSamplingBean bsb = new BailSamplingBean();
        bsb.setReverse(-100, bait_sample_length + 100, true);
        assertThat(bsb.fromStartForward, is(StartForward));
        assertThat(bsb.fromStartReverse, is(StartReverse));
        assertThat(bsb.fromEndForward, is(EndForward));
        assertThat(bsb.fromEndReverse, is(EndReverse));
    }
}
