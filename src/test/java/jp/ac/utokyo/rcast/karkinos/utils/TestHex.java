package jp.ac.utokyo.rcast.karkinos.utils;

import java.math.BigInteger;

import org.apache.commons.codec.DecoderException;
import org.apache.commons.codec.binary.Hex;

public class TestHex {

	public static void main(String[] args) {
		
		//
		String s = "0x0500000a0005000002000100";
		try {
			
			BigInteger bi = new BigInteger(s.substring(2),16);
			byte[] ba = bi.toByteArray();
			for(byte b:ba){
				
				System.out.println(b);
				
			}
			
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}

}
