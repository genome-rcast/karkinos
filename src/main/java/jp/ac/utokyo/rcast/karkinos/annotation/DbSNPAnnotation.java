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
package jp.ac.utokyo.rcast.karkinos.annotation;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.math.BigInteger;
import java.nio.channels.FileChannel;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.codec.binary.Hex;

import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;

public class DbSNPAnnotation implements java.io.Serializable {

	public static final int MODEdbSNP = 0;
	public static final int MODE1000g = 1;
	public static final int MODEcosmic = 2;
	public static final int MODEexonSNP = 3;
	public static final int MODEdbSNPVCF = 4;

	String file = null;
	String g1000 = null;
	String cosmic = null;
	String exonSNP = null;

	float g1000thres = 0.03f;
	float exonThres = 0.03f;

	Map<String, Long> index = null;
	Map<String, Long> indexg1000 = null;
	Map<String, Long> indexcosmic = null;
	Map<String, Long> indexexonSNP = null;

	public DbSNPAnnotation(String _file, String g1000, float g1000thres,
			String cosmic, String exonSNP) {

		//

		if (IsNotEmpty(exonSNP)) {
			try {
				readsIndex(exonSNP, MODEexonSNP);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			this.exonSNP = exonSNP;
		}

		if (IsNotEmpty(g1000)) {
			try {
				readsIndex(g1000, MODE1000g);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			this.g1000 = g1000;
			this.g1000thres = g1000thres;
		}

		if (IsNotEmpty(_file)) {
			
			int mode = MODEdbSNP;
			if(_file.endsWith("vcf")){
				mode = MODEdbSNPVCF;
			}
			try {
				readsIndex(_file,mode);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			file = _file;
		}

		if (IsNotEmpty(cosmic)) {
			try {
				readsIndex(cosmic, MODEcosmic);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			this.cosmic = cosmic;
		}

	}

	private boolean IsNotEmpty(String s) {
		
		return (s!=null) && (s.trim().length() > 2);

	}

	public void readsIndex(String file, int mode) throws IOException {
		// binbitsize fix 23 do not chage for constancy

		boolean vcf = file.endsWith("vcf");
		if (mode == MODEdbSNP) {
			index = new LinkedHashMap<String, Long>();
		} else if (mode == MODEcosmic) {
			indexcosmic = new LinkedHashMap<String, Long>();
		} else if (mode == MODE1000g) {
			indexg1000 = new LinkedHashMap<String, Long>();
		} else {
			indexexonSNP = new LinkedHashMap<String, Long>();
		}
		String indexf = file + ".karkinos.index";

		File f = new File(indexf);
		if (f.exists() && newer(f, new File(file))) {
			// load exsisting file
			BufferedReader br = new BufferedReader(new InputStreamReader(
					new FileInputStream(f)));
			for (;;) {
				String line = br.readLine();
				if (line == null)
					break;

				String[] sa = line.split("\t");
				String key = sa[0];
				long filepos = Long.parseLong(sa[1]);
				if (mode == MODEdbSNP) {
					index.put(key, filepos);
				} else if (mode == MODEcosmic) {
					indexcosmic.put(key, filepos);
				} else if (mode == MODE1000g) {
					indexg1000.put(key, filepos);
				} else {
					indexexonSNP.put(key, filepos);
				}
			}
			br.close();
		} else {

			// FileInputStream fis = new FileInputStream(file);
			// FileChannel fc = fis.getChannel();
			// BufferedReader br = new BufferedReader(new
			// InputStreamReader(fis));

			RandomAccessFile raf = new RandomAccessFile(file, "r");

			int totalcnt = 0;

			String key = "";
			// long filepos =0;
			for (;;) {

				long filepos = raf.getFilePointer();
				String line = raf.readLine();
				if (line == null)
					break;

				// filepos = filepos + line.length()+2;
				// System.out.println(line+filepos);

				if (line.startsWith("#")) {
					continue;
				}

				totalcnt++;
				String[] sa = line.split("\t");
				if (sa.length < 2) {
					continue;
				}
				String chr = null;
				int pos = 0;

				if (mode == MODEdbSNP) {
					chr = sa[1];
					pos = Integer.parseInt(sa[3]);
					if (vcf) {
						chr = sa[0];
						pos = Integer.parseInt(sa[1]);
					}
//					if(sa[10].equals("cDNA")){
//						continue;
//					}
					
				} else {
					chr = sa[0];
					pos = Integer.parseInt(sa[1]);
				}
				if (chr.contains("chr"))
					chr = chr.replace("chr", "");
				String _key = chr + "-" + getBinidx(pos);
				if (!key.equals(_key)) {

					// System.out.println(line+"\t"+filepos);
					if (key.equals("")) {
						_key = chr + "-" + 0;
					}
					if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
						index.put(_key, filepos);
					} else if (mode == MODEcosmic) {
						indexcosmic.put(_key, filepos);
					} else if (mode == MODE1000g) {
						indexg1000.put(_key, filepos);
					} else {
						indexexonSNP.put(_key, filepos);
					}
				}
				key = _key;

			}
			raf.close();

			// write index file
			// if(f.canWrite()){
			// can write does not work in some environments try catch instead

			try {

				BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
						new FileOutputStream(f)));
				Map<String, Long> map = null;
				if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
					map = index;
				} else if (mode == MODEcosmic) {
					map = indexcosmic;
				} else if (mode == MODE1000g) {
					map = indexg1000;
				} else {
					map = indexexonSNP;
				}
				for (Entry<String, Long> es : map.entrySet()) {
					bw.write(es.getKey() + "\t" + es.getValue() + "\n");
				}
				bw.close();
			} catch (IOException ex) {
			}
			// }

		}

	}

	private boolean newer(File f, File file2) {
		return f.lastModified() > file2.lastModified();
	}

	// create index for dbSNP by bin
	public static final int binbitsize = 21;

	private int getBinidx(int pos) {
		return pos >> binbitsize;
	}

	String presentindex = "";
	Map<Integer, DbSNPBean> data = new HashMap<Integer, DbSNPBean>();
	Map<Integer, DbSNPBean> g1000data = new HashMap<Integer, DbSNPBean>();
	Map<Integer, DbSNPBean> cosmicdata = new HashMap<Integer, DbSNPBean>();
	Map<Integer, DbSNPBean> exonSNPdata = new HashMap<Integer, DbSNPBean>();

	public DbSNPBean lookup(String chr, int pos) {

		if (chr.contains("chr")) {
			chr = chr.replace("chr", "");
		}
		String binid = chr + "-" + getBinidx(pos);
		if (!binid.equals(presentindex)) {
			System.out.println("chr" + binid + " " + pos);
			try {
				boolean success = false;
				if (file != null) {
					success = loadData(binid, 0);
					if (success == false)
						return null;
				}				
				if (g1000 != null) {
					success = loadData(binid, 1);
					if (success == false)
						return null;
				}
				if (exonSNP != null) {
					success = loadData(binid, 3);
					if (success == false)
						return null;
				}
				if (cosmic != null) {
					success = loadData(binid, 2);
					if (success == false)
						return null;
				}

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
		}
		DbSNPBean dbSNPbean = data.get(pos);
		if (dbSNPbean != null) {
			dbSNPbean.setMode(MODEdbSNP);
		}
		if (g1000 != null) {
			DbSNPBean dbSNPg1000 = g1000data.get(pos);
			if (dbSNPg1000 != null) {

				if (dbSNPbean == null) {
					dbSNPbean = dbSNPg1000;
					dbSNPbean.setMode(MODE1000g);
				} else {
					if (dbSNPg1000.isValid()) {
						dbSNPbean.setValid(true);
					}
				}

			}

		}
		if (exonSNP != null) {

			DbSNPBean dbSNPexon = exonSNPdata.get(pos);
			if (dbSNPexon != null) {

				if (dbSNPbean == null) {
					dbSNPbean = dbSNPexon;
					dbSNPbean.setMode(MODEexonSNP);
				} else {
					if (dbSNPexon.isValid()) {
						dbSNPbean.setValid(true);
					}
				}

			}

		}
		DbSNPBean cbean = cosmicdata.get(pos);
		if (cbean != null) {
			//
			if (dbSNPbean == null) {
				cbean.setCosmic(true);
				cbean.setCosmicvalid(cbean.valid);
				cbean.setCosmiccount(cbean.cnt);
				return cbean;
			} else {
				dbSNPbean.setCosmic(true);
				dbSNPbean.setCosmicvalid(cbean.valid);
				dbSNPbean.setCosmiccount(cbean.cnt);
			}
		}
		return dbSNPbean;

	}

	private boolean loadData(String binid, int mode) throws IOException {

		Long fp = 0l;
		if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
			data = null;
			data = new HashMap<Integer, DbSNPBean>();
			fp = index.get(binid);

		} else if (mode == MODEcosmic) {
			cosmicdata = null;
			cosmicdata = new HashMap<Integer, DbSNPBean>();
			fp = indexcosmic.get(binid);

		} else if (mode == MODE1000g) {
			g1000data = null;
			g1000data = new HashMap<Integer, DbSNPBean>();
			fp = indexg1000.get(binid);

		} else {

			exonSNPdata = null;
			exonSNPdata = new HashMap<Integer, DbSNPBean>();
			fp = indexexonSNP.get(binid);

		}

		if (fp == null)
			return false;
		FileInputStream fis = null;
		BufferedReader br = null;
		try {
			if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
				fis = new FileInputStream(file);
			} else if (mode == MODEcosmic) {
				fis = new FileInputStream(cosmic);
			} else if (mode == MODE1000g) {
				fis = new FileInputStream(g1000);
			} else {
				fis = new FileInputStream(exonSNP);
			}
			fis.skip(fp);
			br = new BufferedReader(new InputStreamReader(fis));
			int totalcnt = 0;
			long init = fis.available();
			int loadtotal = 0;
			int cntignore = 0;
			for (;;) {
				long filepos = init - fis.available();
				String line = br.readLine();
				if (line == null)
					break;

				totalcnt++;
				String[] sa = line.split("\t");
				if (sa.length < 2) {
					continue;
				}

				String chr = "";
				int pos = 0;
				float freq = 0;
				boolean validated = false;
				String vs = "";

				try {
					chr = null;
					pos = 0;

					if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
						chr = sa[1];
						pos = Integer.parseInt(sa[3]);
						boolean vcf = file.endsWith("vcf");
						if (vcf) {
							chr = sa[0];
							pos = Integer.parseInt(sa[1]);
						}
						try {
							if(mode == MODEdbSNP){
							 vs = sa[12].trim();
							 validated = !(vs.contains("unknown") || vs
									.contains("by-cluster"));
							}
							
							byte[] ba = null;
							if(mode == MODEdbSNPVCF){
								
								String info = sa[5];
								String[] ary = info.split(";");
								String vp ="";
								for(String ss:ary){
									if(ss.contains("VP")){
										vp = ss.replace("VP=","");
										//
										BigInteger bi = new BigInteger(vp.substring(2),16);
										ba = bi.toByteArray();
									}
								
									
								}
								//
								if(ba!=null){
									
									validated = ((ba[6] & 0x4) != 0);
									boolean assemblyproblem  = ((ba[5] & 0x4) != 0) || ((ba[5] & 0x8) != 0);
									if(assemblyproblem){
										validated = false;
									}
									
								}							
															
								
							}
							
							
						} catch (Exception ex) {

						}
					} else if (mode == MODEexonSNP) {

						String id = sa[2];
						if (id.length() > 2) {
							// do not register dbSNP pos
							continue;
						}

						chr = sa[0];
						pos = Integer.parseInt(sa[1]);
						try {
							freq = Float.parseFloat(sa[6]);
//							if (freq < 0.2) {
//								continue;
//							} else if (freq >= exonThres) {
//								validated = true;
//							}
							if (freq >= exonThres) {
								validated = true;
							}
						} catch (Exception ex) {

						}
					} else if (mode == MODEcosmic) {

						chr = sa[0];
						pos = Integer.parseInt(sa[1]);
						validated = false;
						try {

							if(sa.length>8){
							 validated = sa[8].equalsIgnoreCase("true");
							}
							
						} catch (Exception ex) {

						}

					} else {
						
						chr = sa[0];
						pos = Integer.parseInt(sa[1]);
						try {
							freq = Float.parseFloat(sa[4]);
						} catch (Exception ex) {
							
							String info = sa[5];
							String[] ary = info.split(";");
							for(String ss:ary){
								if(ss.contains("AF")){
									String af = ss.replace("AF=","");
									//
									try {
										freq = Float.parseFloat(af);
									} catch (Exception e) {
										
									}
									
								}
							
								
							}
							
						}

					}
					// if((mode ==
					// MODEcosmic)&&(binid.equals("1-4")||binid.equals("1-3"))){
					// System.out.println("check");
					// }

					if (chr.contains("chr")) {
						chr = chr.replace("chr", "");
					}

					if (chr.equals("X") || chr.equals("Y") || chr.equals("M")) {

					} else {

						try {
							int n = Integer.parseInt(chr);
							if (n > 30) {
								continue;
							}
						} catch (Exception ex) {
							continue;
						}

					}

				} catch (Exception ex) {
					continue;
				}
				String key = chr + "-" + getBinidx(pos);
				if (!key.equals(binid)) {
					cntignore++;

					if (cntignore < 5) {
						continue;
					} else {
						break;
					}
				}
				DbSNPBean dbSNP = new DbSNPBean();
				dbSNP.setData(sa);
				dbSNP.setFreq(freq);
				dbSNP.setVaridationStr(vs);

				if (mode == MODEdbSNP || mode == MODEdbSNPVCF) {
					dbSNP.setValid(validated);
					data.put(pos, dbSNP);
				} else if (mode == MODEcosmic) {

					if (cosmicdata.containsKey(pos)) {
						DbSNPBean dbSNP2 = cosmicdata.get(pos);
						dbSNP2.inc();
						if (validated) {
							dbSNP2.setValid(true);
						}
					} else {
						dbSNP.setValid(validated);
						cosmicdata.put(pos, dbSNP);
					}

				} else if (mode == MODE1000g) {
					if (freq != 0) {
						boolean valid = freq >= g1000thres;
						dbSNP.setValid(valid);
					}
					g1000data.put(pos, dbSNP);
				} else {
					dbSNP.setValid(validated);
					exonSNPdata.put(pos, dbSNP);
				}
				loadtotal++;

			}
			// System.out.println(MODEdbSNP+" "+binid+" load="+loadtotal+" "+cntignore);
			presentindex = binid;
		} finally {
			if (br != null) {
				br.close();
				fis.close();
			}
		}
		return true;

	}

}
