/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.rcast.karkinos.utils;

import java.io.IOException;
import java.io.RandomAccessFile;

public class LittleEndian {

	RandomAccessFile raf;
	byte buf[] = null;
	public LittleEndian(RandomAccessFile _raf){
		
		raf = _raf;
		buf = new byte[4];
		
	}
	
	public int readInt() throws IOException{
		
		int i = 0;
		raf.read(buf);
		i = (buf[0] & 0xff)
			+((buf[1] & 0xff)<<8)
			+((buf[2] & 0xff)<<16)
			+((buf[3] & 0xff)<<24);
		return i;
		
	}
	
}
