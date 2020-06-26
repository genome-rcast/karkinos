package jp.ac.utokyo.rcast.karkinos.exec;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;

import org.junit.jupiter.api.Test;

public class CopyNumberIntervalTest {
  @Test
  void chrTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
    };
    assertThat(cni.getChr(), is("chr1"));
  }

  @Test
  void aafbafTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1", 0.1f, 0.2f) {
      private static final long serialVersionUID = 1L;
    };
    assertThat(cni.getAaf(), is(0.1f));
    assertThat(cni.getBaf(), is(0.2f));
  }

  @Test
  void NoSNPTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setNoSNP(4);
    }};
    assertThat(cni.getNoSNP(), is(4));
  }

  @Test
  void allelicTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setAllelic(true);
    }};
    assertThat(cni.isAllelic(), is(true));
  }

  @Test
  void hdeletionTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setHdeletion(false);
    }};
    assertThat(cni.isHdeletion(), is(false));
  }

  @Test
  void startTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setStart(5);
    }};
    assertThat(cni.getStart(), is(5));
  }

  @Test
  void endTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setEnd(8);
    }};
    assertThat(cni.getEnd(), is(8));
  }

  @Test
  void copynumberTest() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setCopynumber(3.2f);
    }};
    assertThat(cni.getCopynumber(), is(3.2f));
  }

  @Test
  void lengthTest1() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setStart(2);
        setEnd(4);
      }
    };
    assertThat(cni.length(), is(3));
  }

  @Test
  void lengthTest2() {
    final CopyNumberInterval cni = new CopyNumberInterval("chr1") {
      private static final long serialVersionUID = 1L;
      {
        setStart(4);
        setEnd(2);
      }
    };
    assertThat(cni.length(), is(0));
  }
}
