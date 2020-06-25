package jp.ac.utokyo.rcast.karkinos.alleliccnv;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.stream.Stream;
import jp.ac.utokyo.rcast.karkinos.exec.CopyNumberInterval;

// TODO: This test class only aims at testing CopyNumberInterval.
public class CNVUtilsTest {
  private static Stream<Arguments> focalampArgs() {
    return Stream.of(
        Arguments.of(true, new CopyNumberInterval("chr1") {{
          setStart(1);
          setEnd(1);
          setCopynumber(6);
        }}),
        Arguments.of(false, new CopyNumberInterval("chr1") {{
          setStart(1);
          setEnd(1000000);
          setCopynumber(6);
        }}),
        Arguments.of(false, new CopyNumberInterval("chr1") {{
          setStart(1);
          setEnd(1);
          setCopynumber(1);
        }}));
  }

  @ParameterizedTest
  @MethodSource("focalampArgs")
  void focalampTest(final boolean expected, final CopyNumberInterval cni)
      throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
    final Method focalamp = CNVUtils.class.getDeclaredMethod("focalamp", CopyNumberInterval.class);
    focalamp.setAccessible(true);
    assertEquals(expected, (boolean)focalamp.invoke(CopyNumberInterval.class, cni));
  }
}
