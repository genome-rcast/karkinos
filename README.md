[![test](https://github.com/genome-rcast/karkinos/workflows/test/badge.svg)](https://github.com/genome-rcast/karkinos/actions)
[![codecov](https://codecov.io/gh/genome-rcast/karkinos/branch/master/graph/badge.svg)](https://codecov.io/gh/genome-rcast/karkinos)

# About

karkinos is tumor genotyper which detects single nucleotide variation (SNV),
integer copy number variation (CNV) and calculates tumor cellularity from tumor-normal paired sequencing data.
Accurate CNV calling is achieved using continuous wavelet analysis and multi-state HMM,
while SNV call is adjusted by tumor cellularity and filtered by heuristic filtering algorithm and Fisher Test.
Also, Noise calls in low depth region are removed using EM algorithm.


# Licence

 Copyright (C) 2014 Hiroki Ueda Rcast, the University of Tokyo

 Licensed under the Apache License, Version 2.0 (the &quot;License&quot;);
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an &quot;AS IS&quot; BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.




# Install

- Download karkinos.jar from web sites.
- Install Java Runtime 1.5 or later.

# Required files

- Bam files for normal reads Alignment.
- Bam files for tumor reads Alignment.
- CaptureTarget Information (bed format)
- dbSNP (bin,chr,start,end)
- cosmic(vcf)
- karkinos.property

dbSNP file can be download from
https://usegalaxy.org/library_common/ldda_info?library_id=e660f0c750f4b341&show_deleted=False&cntrller=library&folder_id=b89fca8d2667b9c7&use_panels=False&id=c90f779d48eddc3d

where format is following
1 = bin (for indexing)
2 = chromosome
3 = start (0 based)
4 = end (1 based)
5 = rs#
6 = score
7 = strand
8 = ref allele from NCBI
9 = ref allele from UCSC
10 = observed alleles
11 = molType
12 = class
13 = valid
14 = avHet
15 = avHetSE
16 = func
17 = location type
18 = weight
19 = exceptions
20 = submitterCount
21 = submitters
22 = allele frequency count
23 = alleles
24 = alleleNs
25 = alleleFreqs
26 = bitfields

example line

```
585	chr1	10468	10469	rs117577454	0	+	C	C	C/G	genomic	single	by-1000genomes	0	0	unknown	exact	1		1	1000GENOMES,	2	G,C,	18.000000,102.000000,	0.150000,0.850000,
```

# Run

command

- analysis - pileup and analysis SNV,CNV,Tumor purity
- analysis - pileup only if chr is disgnated
- reanalysis - analysis SNV,CNV,Tumor purity from piled up files.

For single computational node use analysis command to excute all genomic range.

If parallel computational nodes are avilavle, use analysis command to each chromosome and use reanalysis command to maerge the result.

Example shell script is on the web sites.

```
karkinos.jar analysis

usage: karkinos.jar analysis -n <arg> -t <arg> -r <arg> -snp <arg> -ct
       <arg> -o <arg> -id <arg> [-prop <arg>] [-mp <arg>] [-g1000 <arg>]
       [-cosmic <arg>] [-g1000freq <arg>] [-chr <arg>] [-rs <arg>] [-rg
       <arg>] [-exonSNP <arg>] [-nopdf]
 -n,--normalBam <arg>                normal bam file
 -t,--tumorBam <arg>                 tumor bam file
 -r,--reference <arg>                2 bit genome reference file
 -snp,--dbSNP <arg>                  dbSNP list from annover
                                     sites,(bin,chr,start,end)
 -ct,--captureTarget <arg>           Capture target regions(bed format)
 -o,--outdir <arg>                   output directory
 -id,--uniqueid <arg>                unique id for this sample
 -prop,--property <arg>              path to property file( otherwise
                                     default val)
 -mp,--mappability <arg>             optional,mappability from ucsc (bw,
                                     big wig format)
 -g1000,--1000genome <arg>           optional,1000 genome list from
                                     annover
                                     sites,(chr,pos,ref,alt,freq,id)
 -cosmic,--cosmicSNV <arg>           cosmic snv vcf format
 -g1000freq,--1000genomefreq <arg>   optional,1000 genome frequency
                                     threshold to use
 -chr,--chrom <arg>                  chromosome to analyze
 -rs,--readsStats <arg>              optional,reads stats
                                     files(normal,tumor)
 -rg,--refFlatGenes <arg>            optional,gene reference for depth
                                     stats
 -exonSNP,--exonSNP <arg>            additional Exon SNP
 -nopdf,--nopdf                      no graphic summary pdf output
```

```
usage: karkinos.jar reanalysis -md <arg> -t <arg> -r <arg> -snp <arg> -ct
       <arg> -o <arg> -id <arg> [-prop <arg>] [-mp <arg>] [-g1000 <arg>]
       [-cosmic <arg>] [-g1000freq <arg>] [-rs <arg>] [-rg <arg>] [-tc
       <arg>] [-exonSNP <arg>] [-nd <arg>] [-nopdf]
 -md,--middledate <arg>              middle data object
 -t,--tumorBam <arg>                 tumor bam file
 -r,--reference <arg>                2 bit genome reference file
 -snp,--dbSNP <arg>                  dbSNP list from annover
                                     sites,(bin,chr,start,end)
 -ct,--captureTarget <arg>           Capture target regions(bed format)
 -o,--outdir <arg>                   output directory
 -id,--uniqueid <arg>                unique id for this sample
 -prop,--property <arg>              path to property file( otherwise
                                     default val)
 -mp,--mappability <arg>             optional,mappability from ucsc (bw,
                                     big wig format)
 -g1000,--1000genome <arg>           optional,1000 genome list from
                                     annover
                                     sites,(chr,pos,ref,alt,freq,id)
 -cosmic,--cosmicSNV <arg>           cosmic snv vcf format
 -g1000freq,--1000genomefreq <arg>   optional,1000 genome frequency
                                     threshold to use
 -rs,--readsStats <arg>              reads stats files(normal,tumor)
 -rg,--refFlatGenes <arg>            optional,gene reference for depth
                                     stats
 -tc,--tumorContents <arg>           fixed tumor contents (and skip tc
                                     calculation)
 -exonSNP,--exonSNP <arg>            additional Exon SNP
 -nd,--normaldepth <arg>             use averaged normal depth if
                                     available
 -nopdf,--nopdf                      no graphic summary pdf output
```