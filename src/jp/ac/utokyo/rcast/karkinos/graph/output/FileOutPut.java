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
package jp.ac.utokyo.rcast.karkinos.graph.output;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import au.com.bytecode.opencsv.CSVReader;
import jp.ac.utokyo.karkinos.noisefilter.NoiseAnalysis;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.AllelicCNV;
import jp.ac.utokyo.rcast.karkinos.alleliccnv.SNVHolderPlusACnv;
import jp.ac.utokyo.rcast.karkinos.annotation.DbSNPAnnotation;
import jp.ac.utokyo.rcast.karkinos.exec.CapInterval;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.exec.PileUPResult;
import jp.ac.utokyo.rcast.karkinos.exec.SNVHolder;
import jp.ac.utokyo.rcast.karkinos.filter.DefinedSites;
import jp.ac.utokyo.rcast.karkinos.readssummary.GeneExons;
import jp.ac.utokyo.rcast.karkinos.readssummary.Interval;
import jp.ac.utokyo.rcast.karkinos.readssummary.ReadsSummary;
import jp.ac.utokyo.rcast.karkinos.utils.DataHolder;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import jp.ac.utokyo.rcast.karkinos.wavelet.WaveletIF;

public class FileOutPut {

	public static void outPutSNVDataVCF(String outpath, DataSet dataset,
			TwoBitGenomeReader tgr, NoiseAnalysis na) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));
			bw.append("##fileformat=VCFv4.1" + "\n");

			bw.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
					+ "\n");
			bw.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency adjusted by tumor ratio\">"
					+ "\n");
			bw.append("##INFO=<ID=AFO,Number=A,Type=Float,Description=\"Allele Frequency original\">"
					+ "\n");
			bw.append("##INFO=<ID=CN,Number=A,Type=Float,Description=\"Copy Number\">"
					+ "\n");
			bw.append("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">"
					+ "\n");
			bw.append("##INFO=<ID=cosmic,Number=0,Type=Flag,Description=\"cosmic DB membership\">"
					+ "\n");

			bw.append("##INFO=<ID=S=0,Type=Float,Description=\"Sequence Entropy\">"
					+ "\n");
			bw.append("##INFO=<ID=MP=0,Type=Float,Description=\"Mappability\">"
					+ "\n");
			bw.append("##INFO=<ID=pD=0,Type=Float,Description=\"Fisher test pval for reads direction support reads vs ref reads\">"
					+ "\n");

			bw.append("##INFO=<ID=pval,Number=0,Type=Float,Description=\"Fisher test pval\">"
					+ "\n");
			bw.append("##INFO=<ID=OD,Number=0,Type=Flag,Description=\"support by only one directinal reads\">"
					+ "\n");
			bw.append("##INFO=<ID=RP,Number=0,Type=Flag,Description=\"masked repeat\">"
					+ "\n");
			bw.append("##INFO=<ID=LowRefFilter,Number=0,Type=Flag,Description=\"low odds log ratio normal\">"
					+ "\n");
			bw.append("##INFO=<ID=LowTumorFilter,Number=0,Type=Flag,Description=\"low odds log ratio tumor\">"
					+ "\n");

			bw.append("##FILTER=<ID=qf,Description=\"Quality below threshold\">"
					+ "\n");
			bw.append("##FILTER=<ID=bf,Description=\"Bayesian filterling\">"
					+ "\n");
			bw.append("##FILTER=<ID=snp,Description=\"dbSNP snp\">" + "\n");
			bw.append("##FILTER=<ID=ma,Description=\"Low mappability \">"
					+ "\n");
			bw.append("##FILTER=<ID=mq,Description=\"Low mapping quality \">"
					+ "\n");
			bw.append("##FILTER=<ID=entropy,Description=\"Low complexty\">"
					+ "\n");
			bw.append("##FILTER=<ID=srd,Description=\"Support reads have only one direction\">"
					+ "\n");
			bw.append("##FILTER=<ID=srm,Description=\"Only support reads have another mismatch\">"
					+ "\n");
			bw.append("##FILTER=<ID=sre,Description=\"mutation at reads ends\">"
					+ "\n");
			bw.append("##FILTER=<ID=scc,Description=\"mutation at same cycle\">"
					+ "\n");
			bw.append("##FILTER=<ID=mmt,Description=\"too many mismach in support reads\">"
					+ "\n");

			bw.append("##INFO=<ID=ND,Number=1,Type=Integer,Description=\"Total Depth Normal\">"
					+ "\n");
			bw.append("##INFO=<ID=NR,Number=1,Type=Integer,Description=\"ratio Normal\">"
					+ "\n");

			// bw.append("##FILTER=<ID=logn,Description=\"logn, odds log ratio normal\">"+
			// "\n");
			// bw.append("##FILTER=<ID=logt,Description=\"logt, odds log ratio tumor\">"+
			// "\n");
			bw.append("##FILTER=<ID=logn,Description=\"logn, odds log ratio \">"
					+ "\n");
			bw.append("##FILTER=<ID=logt,Description=\"logt, odds log ratio tumor adjusted by tumor ratio and Copy number\">"
					+ "\n");
			bw.append("##FILTER=<ID=fisher,Description=\"failed by fisher test\">"
					+ "\n");

			bw.append("##FILTER=<ID=near_indel,Description=\"near_indel\">"
					+ "\n");
			bw.append("##FILTER=<ID=low_adj_ratio ,Description=\"Low_tumor_adjustedRatio\">"
					+ "\n");
			bw.append("##FILTER=<ID=low_adj_reads ,Description=\"Low_tumor_adjustedRead\">"
					+ "\n");

			bw.append("##FILTER=<ID=high_lqr ,Description=\"Too many LowQualReads\">"
					+ "\n");
			bw.append("##FILTER=<ID=high_TCn ,Description=\"HighnormalRatio at TC adjusted cand site\">"
					+ "\n");
			bw.append("##FILTER=<ID=lowTCreads ,Description=\"Low_tumor_adjustedRead\">"
					+ "\n");
			bw.append("##FILTER=<ID=lowTCreads2,Description=\"Low_tumor_adjustedRead\">"
					+ "\n");

			bw.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
					+ "\n");
			bw.append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
					+ "\n");
			bw.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
					+ "\n");

			bw.append("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FILTER2	Normal	Tumor	seqbefore	seqafter"
					+ "\n");

			float tumorRratio = dataset.getTumorRatio();
			for (SNVHolder holder : dataset.getSnvlist()) {
				String str = FormatHelper.getVCFLine(holder, tgr, tumorRratio,
						na);
				if (str != null) {
					bw.write(str + "\n");
				}
			}

			bw.close();
		} catch (IOException ex) {
		}

	}

	public static void outPutCNVData(String outpath, DataSet dataset) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));
			for (List<WaveletIF> list : dataset.getCapInterval()) {
				for (WaveletIF wif : list) {

					CapInterval ci = (CapInterval) wif;
					bw.write(ci.getInfoStr() + "\n");

				}
			}
			bw.close();
		} catch (IOException ex) {
		}

	}

	public static void lowcovBed(String outpathnormal, String outpathtumor,
			ReadsSummary readsSummary, GeneExons ge) {

		//
		List<Interval> lowcoverageListN = readsSummary.getNormalDepth()
				.getLowcoverageList();
		List<Interval> lowcoverageListT = readsSummary.getTumorDepth()
				.getLowcoverageList();

		writeLowDepthBed(outpathnormal, lowcoverageListN,ge);
		writeLowDepthBed(outpathtumor, lowcoverageListT,ge);
	}

	private static void writeLowDepthBed(String filepath,
			List<Interval> lowcoverageList, GeneExons ge) {
		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(filepath)));
			bw.append("#chr\t	start\t	end\t	name\t	depth\t geneSymbol\t  refseqid\t exonIdx \n");
			for (Interval iv : lowcoverageList) {
				bw.append(iv.getChr() + "\t" + iv.getStart() + "\t"
						+ iv.getEnd() + "\t" + iv.getChr() + ":"
						+ iv.getStart() + "-" + iv.getEnd() + "\t"
						+ iv.getDepth());

				if (ge != null) {

					Interval iv2 = ge.getGeneIntervalExon(iv.getChr(),
							iv.getStart());
					if (iv2 == null) {
						iv2 = ge.getGeneIntervalExon(iv.getChr(), iv.getEnd());
					}
					//
					if (iv2 == null) {
						bw.write("\t \t \t");
					} else {
						bw.write("\t" + iv2.getGeneSymbol());
						bw.write("\t" + iv2.getRefseqid());
						bw.write("\t" + iv2.getExonidx());
					}

				} else {
					bw.write("\t \t");
				}

				bw.write("\n");

			}
			bw.close();
		} catch (IOException ex) {
		}
	}

	public static void outPutSNVDataForAnnover(String outpath, DataSet dataset,
			TwoBitGenomeReader tgr) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			float tumorRratio = dataset.getTumorRatio();
			for (SNVHolder holder : dataset.getSnvlist()) {
				String str = FormatHelper.getAnnoverInputLine(holder, tgr,
						tumorRratio);
				if (str != null) {
					bw.write(str + "\n");
				}
			}

			bw.close();
		} catch (IOException ex) {
		}

	}

	// chr
	// pos
	// id
	// ref
	// alt
	// qual
	// refcnt;
	// altcnt;
	// freq;
	// alt
	// qual
	// refcnt;
	// altcnt;
	// freq;
	// info
	// tn ratio

	public static void sites(String outpath, DataSet dataset,
			TwoBitGenomeReader tgr, GeneExons ge, String sites) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			DefinedSites ds = new DefinedSites(sites);

			bw.append("#CHROM	POS	ID	REF	ALT_N	QUAL_N	REFCNT_N	ALTCNT_N 	FEEQ_N	"
					+ "ALT_T	QUAL_T	REFCNT_T	ALTCNT_T 	FEEQ_T	INFO	TNRATIO"
					+ "\n");

			// float tumorRratio = dataset.getTumorRatio();
			for (SNVHolder holder : dataset.getSnvlist()) {

				String chr = holder.getChr();
				int pos = holder.getPos();
				if (ds.contains(chr, pos)) {

					String[] stra = FormatHelper.getLine(holder, tgr, ge);
					if (stra != null) {

						for (String s : stra) {

							bw.write(s);
							bw.write("\t");
						}
						bw.write("\n");
					}

				}
			}
			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public static void allDiff(String outpath, DataSet dataset,
			TwoBitGenomeReader tgr, GeneExons ge) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));
			bw.append("##fileformat=VCFv4.1" + "\n");

			bw.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
					+ "\n");
			bw.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
					+ "\n");
			bw.append("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">"
					+ "\n");
			bw.append("##INFO=<ID=geneID,Number=0,Type=Flag,Description=\"refseq gene id\">"
					+ "\n");
			bw.append("##INFO=<ID=onCDS,Number=0,Type=Flag,Description=\"on cds\">"
					+ "\n");
			bw.append("##INFO=<ID=1000freq,Number=0,Type=Flag,Description=\"allele frequency of 1000 genome\">"
					+ "\n");
			bw.append("#CHROM	POS	ID	REF	ALT	QUAL		" + "\n");

			// float tumorRratio = dataset.getTumorRatio();
			for (SNVHolder holder : dataset.getSnvlist()) {

				String str = FormatHelper.getVCFLineAllDiff(holder, tgr, ge);
				if (str != null) {
					bw.write(str + "\n");
				}
			}
			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public static void outputSNP(String outpath, DataSet dataset,
			TwoBitGenomeReader tgr, GeneExons ge) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));
			bw.append("##fileformat=VCFv4.1" + "\n");

			bw.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
					+ "\n");
			bw.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
					+ "\n");
			bw.append("##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">"
					+ "\n");
			bw.append("##INFO=<ID=geneID,Number=0,Type=Flag,Description=\"refseq gene id\">"
					+ "\n");
			bw.append("##INFO=<ID=onCDS,Number=0,Type=Flag,Description=\"on cds\">"
					+ "\n");
			bw.append("##INFO=<ID=1000freq,Number=0,Type=Flag,Description=\"allele frequency of 1000 genome\">"
					+ "\n");
			bw.append("#CHROM	POS	ID	REF	ALT	QUAL		" + "\n");

			// float tumorRratio = dataset.getTumorRatio();
			for (SNVHolder holder : dataset.getSnvlist()) {
				String str = FormatHelper.getVCFLine4SNP(holder, tgr, ge);
				if (str != null) {
					bw.write(str + "\n");
				}
			}
			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public static void outputDepth(DataSet dataset, String outpath,
			float purity, GeneExons ge) {

		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			bw.write("#chr" + "\t" + "start" + "\t" + "end" + "\t"
					+ "normaldepth" + "\t" + "tumordepth" + "\t" + "t/nRatio"
					+ "\t" + "t/nRatioAdj" + "\t" + "smoothedRatio" + "\t"
					+ "adjustedRatio" + "\t" + "AalleleCopyNumber" + "\t"
					+ "BalleleCopyNumber" + "\t" + "copyNumber" + "\t" + "gene"
					+ "\t" + "id" + "\t" +"exonidx"+ "\n");

			float factor = 1 / purity;
			for (List<WaveletIF> list : dataset.getCapInterval()) {

				for (WaveletIF wi : list) {

					CapInterval ci = (CapInterval) wi;
					bw.write(ci.getChr() + "\t" + ci.getStart() + "\t"
							+ ci.getEnd() + "\t"
							+ format(ci.getCNVInfo().getNormaldepthAdj())
							+ "\t" + format(ci.getCNVInfo().getTumordepthAdj())
							+ "\t"
							+ format(ci.getCNVInfo().getOriginalTnratio())
							+ "\t" + format(ci.getCNVInfo().getTnratio())
							+ "\t" + format(ci.getCNVInfo().getDenoise())
							+ "\t"
							+ format(ci.getCNVInfo().getDenoise() * factor)
							+ "\t" + format(ci.getAafreq()) + "\t"
							+ format(ci.getBafreq()) + "\t"
							+ format(ci.getHMMValue() * 2));

					if (ge != null) {

						Interval iv = ge.getGeneIntervalExon(ci.getChr(),
								ci.getStart());
						if (iv == null) {
							iv = ge.getGeneIntervalExon(ci.getChr(), ci.getEnd());
						}
						//
						if (iv == null) {
							bw.write("\t \t");
						} else {
							bw.write("\t" + iv.getGeneSymbol());
							bw.write("\t" + iv.getRefseqid());
							bw.write("\t" + iv.getExonidx());
						}

					} else {
						bw.write("\t \t \t");
					}

					bw.write("\n");

				}
			}

			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	private static String format(double num) {
		NumberFormat nf = NumberFormat.getNumberInstance();
		return nf.format(num);
	}

	public static void outputAlleleDepth(AllelicCNV alCNV, String outpath,
			float purity) {
		try {

			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			bw.write("#chr" + "\t" + "pos" + "\t" + "id" + "\t"
					+ "highAlleleRow" + "\t" + "lowAlleleRow" + "\t"
					+ "highAlleleAdj" + "\t" + "lowAlleleAdj" + "\t"
					+ "highSmooth" + "\t" + "lowSmooth" + "\t"
					+ "highSmoothAdj" + "\t" + "lowSmoothAdj" + "\t"
					+ "highHMM" + "\t" + "lowHMM" + "\n");

			float factor = 1 / purity;
			for (List<SNVHolderPlusACnv> plist : alCNV.getList()) {

				for (SNVHolderPlusACnv sc : plist) {

					String id = ".";
					if (sc.getSnv().getDbSNPbean() != null) {
						// id
						id = sc.getSnv().getDbSNPbean().getInfo();

					}

					bw.write(sc.getSnv().getChr() + "\t" + sc.getSnv().getPos()
							+ "\t" + id + "\t" + sc.getHighera().getRow()
							+ "\t" + sc.getLowera().getRow() + "\t"
							+ sc.getHighera().getGcadjusted() + "\t"
							+ sc.getLowera().getGcadjusted() + "\t"
							+ sc.getHighera().getWtval() + "\t"
							+ sc.getLowera().getWtval() + "\t"
							+ sc.getHighera().getWtval() * factor + "\t"
							+ sc.getLowera().getWtval() * factor + "\t"
							+ sc.getHighera().getHmmval() + "\t"
							+ sc.getLowera().getHmmval() + "\n");

				}
			}

			bw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public static final String title = "genesymbol \t name \t chrom \t strand \t txStart \t txEnd"
			+ "\t cdsStart \t cdsEnd \t ExonCnt \t starts \t ends "
			+ "\t normalcnt \t tumorcnt \t normalNorm \t tumorNorm \t tnratio \t LRR "
			+ "\t tnsmooth \t cnvhmm \t  AAF \t BAF";

	public static void outputgeneCNV(String outpath, GeneExons ge,
			String refflat) {

		try {
			//
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(outpath)));

			bw.write(title);
			bw.write("\n");

			Map<String, DataHolder> counterForGene = ge.getCounterForGene();

			CSVReader brcvs = new CSVReader(new FileReader(refflat), '\t');

			String[] data = null;
			List<String> l = null;

			while ((data = brcvs.readNext()) != null) {

				String refseqid = data[1];
				l = new ArrayList<String>();
				for (String d : data)
					l.add(d);
				//
				DataHolder df = counterForGene.get(refseqid);
				if (df != null) {
					df.setTotal(ge.getNormaltotal(), ge.getTumortotal());

					l.add(String.valueOf(df.getNormalcnt()));
					l.add(String.valueOf(df.getTumorcnt()));

					l.add(String.valueOf(df.normalAdjust()));
					l.add(String.valueOf(df.tumorAdjust()));

					l.add(String.valueOf(df.tumorRatio()));
					l.add(String.valueOf(df.logr()));

					l.add(String.valueOf(df.averages()[0]));
					l.add(String.valueOf(df.averages()[1]));
					l.add(String.valueOf(df.averages()[2]));
					l.add(String.valueOf(df.averages()[3]));
					l.add(String.valueOf(df.averages()[4]));

				} else {

					l.add(String.valueOf(0));
					l.add(String.valueOf(0));

					l.add(String.valueOf(0));
					l.add(String.valueOf(0));

					l.add(String.valueOf(0));
					l.add(String.valueOf(0));

					l.add(String.valueOf(0));
					l.add(String.valueOf(0));
					l.add(String.valueOf(0));
					l.add(String.valueOf(0));
					l.add(String.valueOf(0));

				}

				boolean first = true;
				for (String s : l) {

					if (first) {
						first = false;

					} else {
						bw.write("\t");
					}
					bw.write(s);

				}
				bw.write("\n");
			}

			bw.close();
			brcvs.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

}
