# Output description for panelLowFreq pipeline

**panelLowFreq** is a bioinformatics best-practice variant calling analysis pipeline used for WES-Seq (whole exome sequencing) or target sequencing. The pipeline focused in variant calling and annotation of candidate low frequency variants.

This document describes the output produced by the pipeline and location of output files.

## Pipeline overview:
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Trimmomatic](#trimming) - adapter and low quality trimming
* [BWA](#bwa) - mapping against reference genome
* [SAMtools](#samtools) - alignment result processing
* [Picard](#picard) - enrichment and alignment metrics
* [GATK](#varscan) - variant calling.
* [KGGSeq](#kggseq) - variant annotation.
* [MultiQC](#multiqc) - quality statistics summary

> Each analysis folder contains a log folder with the log files for each process and each sample.
## Preprocessing
### FastQC
Quality control is performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.
For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Results directory**: ANALYSIS/{ANALYSIS_ID}/QC/fastqc
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_fastqc.html`: html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
   - `{sample_id}/{sample_id}_R[12]_fastqc` : folder with fastqc output in plain text.
   - `{sample_id}/{sample_id}_R[12]_fastqc.zip`: zip compression of above folder.
   - `{sample_id}/{sample_id}.trimmed_R[12]_fastqc.html`: html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC for trimmed reads.
   - `{sample_id}/{sample_id}.trimmed_R[12]_fastqc`: folder with fastqc output in plain text for trimmed reads.
   - `{sample_id}/{sample_id}.trimmed_R[12]_fastqc.zip`: zip compression of above folder for trimmed reads.

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used for removal of adapter contamination and trimming of low quality regions. 
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 15 in a 4 nucleotide window.
-  Read lenght < 70

MultiQC reports the percentage of bases removed by trimming in bar plot showing percentage or reads trimmed in forward and reverse.

**Note**:The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `ANALYSIS/{ANALYSIS_ID}/QC/fastqc` directory.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/QC/trimmomatic
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12].trimmed.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
   - `{sample_id}/{sample_id}_R[12]_unpaired.fastq.gz`: contains high quality reads with only forward or reverse tags surviving.

**NOTE:** This results are not delivered to the researcher by default due to disk space issues. If you are interesested in using them, please contact us and we will add them to your delivery.

## Mapping
### BWA
[BWA](http://bio-bwa.sourceforge.net/), or Burrows-Wheeler Aligner, is designed for mapping low-divergent sequence reads against reference genomes. The result alignment files are further processed with [SAMtools](http://samtools.sourceforge.net/), sam format is converted to bam, sorted and an index *.bai* is generated.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/Alignment
- There is one folder per sample.
- This files can be used in [IGV](https://software.broadinstitute.org/software/igv/) for alignment visualization.
- Files:
   - `{sample_id}/{sample_id}_sorted.bam` : sorted aligned bam file.
   - `{sample_id}/{sample_id}_sorted.bam.bai`: index file for soreted aligned bam.

**NOTE:** This results are not removed due to disk space issues, only last bam processed bam file is retained.. If you are interesested in using them, please contact us and we will try to generate them and add them to your delivery.

## Picard
### MarkDuplicates
The [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) module in the [Picard](https://broadinstitute.github.io/picard/) toolkit differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores, which helps to identify duplicates that arise during sample preparation e.g. library construction using PCR.

The Picard section of the MultiQC report shows a bar plot with the numbers and proportions of primary reads, duplicate reads and unmapped reads.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/Alignment
- There is one folder per sample.
- This files can be used in [IGV](https://software.broadinstitute.org/software/igv/) for alignment visualization.
- Files:
   - `{sample_id}/{sample_id}.woduplicates.bam` : sorted aligned bam file.
   - `{sample_id}/{sample_id}.woduplicates.bam.bai`: index file for soreted aligned bam.

### HsMetrics
Metrics for the analysis of target-capture sequencing experiments are calculated with [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics). The metrics in this class fall broadly into three categories:

- Basic sequencing metrics that are either generated as a baseline against which to evaluate other metrics or because they are used in the calculation of other metrics. This includes things like the genome size, the number of reads, the number of aligned reads etc.
- Metrics that are intended for evaluating the performance of the wet-lab assay that generated the data. This group includes metrics like the number of bases mapping on/off/near baits, %selected, fold 80 base penalty, hs library size and the hs penalty metrics. These metrics are calculated prior to some of the filters are applied (e.g. low mapping quality reads, low base quality bases and bases overlapping in the middle of paired-end reads are all counted).
- Metrics for assessing target coverage as a proxy for how well the data is likely to perform in downstream applications like variant calling. This group includes metrics like mean target coverage, the percentage of bases reaching various coverage levels, and the percentage of bases excluded by various filters. These metrics are computed using the strictest subset of the data, after all filters have been applied.

The Picard section of the MultiQC report shows a table with data as the Bait design efficiency, Bait territory, Fold 80 base penalty, Fold enrichment; near, on and off bait bases; on target bases, Pf (passing filter) reads, Pf unique reads, Pf uq (unique) bases aligned, Pf uq reads aligned, Total reads	Mean bait coverage, Mean target coverage, On bait vs selected, Target territory, Zero cvg targets pct.

**Results directory:** ANALYSIS/{ANALYSIS_ID}/stats/bamstats
- Files:
   - `hsMetrics_all.out` : summary of the some of the most meaningful columns in picard hsmetrics output for all the samples in the project.
   - `{sample_id}_hsMetrics.out`: full picard hsmetrics output per sample.
   - Description of Picard hsMetrics columns in its output can be found in [AnnexIII](#annex-iii)
   
## Variant Calling
### GATK
Best Practice GATK protocol has been used for the variant calling of germinal variations in the family (http://www.broadinstitute.org/gatk/guide/best-practices, november 2015).
This workflow comprises three different steps:

1. **Data preprocessing:**
  a.  **Realignment:** starting from BAM file generated with [BWA](#bwa) and [Picard MarkDuplicates](#markduplicates), realignment around candidate indels is performed in order to improve mapping in complicated zones (low complexity, homopolymers, etc). 
  
   **Results directory**: ANALYSIS/{ANALYSIS_DIR}/variants/variants_gatk/realignment
    - `{sample_id}.realigned.bam` : realigned bam.
    - `{sample_id}.realigned.bam-bai`: index for realigned bam.
    - `{sample_id}.woduplicates.bam-IndelRealigner.intervals`: file with "problematic" regions where realignment was needed.
    
   **NOTE:** This results are not removed due to disk space issues, only last bam processed bam file is retained.. If you are interesested in using them, please contact us and we will try to generate them and add them to your delivery.
  
  b. **Base Recalibrarion**: Next step carries out a phed quality recalibration of bases, using a gold standard set of known SNPS (dbSNP138) using a machine learning approach.
   
   **Results directory**: ANALYSIS/{ANALYSIS_DIR}/variants/variants_gatk/recalibration
    - `{sample_id}.woduplicates.bam-BQSR.pdf` : pdf file with quality graphics before and after recalibration.
    - `{sample_id}.recalibrated.bam` : recalibrated bam.
    - `{sample_id}.recalibrated.bai` : index for recalibrated bams.
    - `{sample_id}.woduplicates.bam-BQSR.csv`: intermediate file.
    - `{sample_id}.woduplicates.bam-recal2_data.grp`: intermediate file.
    
   **NOTE:** BAM files are removed due to disk space issues. If you are interesested in using them, please contact us and we will try to generate them and add them to your delivery.

2. **Variant Calling:**

Se realiza la llamada a variantes de forma conjunta para los tres individuos de la familia. El formato saca el genotipo de cada individuo para todas las posiciones en las que alguno de ellos tenga una variante respecto a referencia. Se utiliza el nuevo módulo de variant calling de GATK, HaplotypeCaller que se vale del cálculo de haplotipos para mejorar la llamada de variantes. Además, es capaz de llamar SNPs, indels y algunas variantes estructurales haciendo un ensamblado de novo. 

Determina si una región es potencialmente variable
Construye un ensamblado deBruigin de la región.
Los paths en el grafo son haplotipos potenciales que tienen que ser evaluados.
Se calcula los likelihoods de los haplotipos dados los datos usando un modelo PairHMM.
Determina si hay alguna vairante ente los haplotipos más probables.
Calcula la distribución de la frecuencia alélica para determinar el contaje de alelos más probables y emite una variante si se da el caso.
Si se emite una variante se calcula el genotipo para cada muestra.

Se analiza con una diferencia respecto a las Best Practice de GATK. Se analizan las 5 muestras al mismo tiempo con HaplotypeCaller (esto cambiará en las siguientes versiones donde se analizarán por separado generando gVCFs en vez de VCFs, ya que en este modo ya se hace el phasing directamente con HaplotypeCaller). Se obtiene por tanto un vcf con todas las variantes detectadas en alguno de los 5 miembros.

HardFiltering: este paso se hace ante la imposibilidad de realizar una recalibración de variantes según el workflow de GATK (VariantRecalibration) que necesita de al menos 10 muestras para que los scores estadísticos sean significativos. Por lo tanto, tenemos que realizar filtros simples para quedarnos con las variantes de mejor calidad y eliminar los posibles falsos positivos según los filtros recomendados en las Best Practices de GATK:

MQ < 40. RMSMappingQuality. This is the Root Mean Square of the mapping quality of the reads across all samples.
DP <5. LowCoverage
QD <2.0. LowQD. This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-reference samples.
FS >60.0. p-value StrandBias. Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls
MQRankSum < -12.5. MappingQualityRankSumTest. This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
ReadPosRankSum < -8.0. VariantReadPosEnd. This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test cannot be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
SOR > 4.0. StrandOddRank. Strandbias.

-	Refinamiento de los genotipos: una vez que se obtienen los genotipos para todas las posiciones variantes de la familia, hay dos herramientas de gatk que se pueden aplicar para refinar los datos obtenidos sabiendo el pedigrí. Estos son:

o	PhaseByTransmission: técnica estadística para determinar cuándo es posible qué alelo proviene del padre y cuál de la madre en el niño. Por consenso se pone primero el alelo de la madre y luego el del padre. Madre 1/0 padre 0/0, niño 1(madre) | 0(padre). Sólo se permite el análisis de trios actualmente, por lo que sólo aparecerán las fases calculadas para padre, madre y probando. No se realizará phasing en los hermanos.
o	ReadBackedPhasing: determina la presencia de haplotipos en cada muestra, no entre ellas. Busca grupos de snps que se encuentran en el mismo cromosoma. Por lo que entiendo haplotypeCaller utiliza este mecanismo para determinar los genotipos pero no lo marca en el vcf. Al correr este Walker en el vcf se marcan los genotipos con la | en vez de con / cuando se ha conseguido determinar un haplotipo.

Los resultados del variant calling se encuentran en ..\ANALYSIS\20180220_TRIO04401 \variant_calling\variants_gatk\variants. Ahí están los ficheros de variantes, los ficheros donde se encuentran las variantes ya filtradas y anotados los haplotipos ..\ANALYSIS\20180220_TRIO04401\variant_calling\variants_gatk\variants\all_sample_varants_gtpos_fil.vcf. Este último fichero será el que hay que anotar for gen, efecto génico, buscar los snps que sigan el modelo de la enfermedad (de novo, heterocigosis compuesta)

Post-Análisis: Anotación y filtrado
Para el post análisis de las variantes obtenidas con GATK se utiliza el software KGGSeq (Li, Gui, Kwan, Bao, & Sham, 2012), una herramienta diseñada para priorizar variantes en el estudio de enfermedades mendelianas. En resumen, se trata de una herramienta de anotación de variantes. Permite incluir información de efecto, gen, tránscrito, enfermedades conocidas relacionadas, artículos de pubmed, etc. Además, permite realizar distintos tipos de filtro por calidad phred, por frecuencia, por common_snps (con frecuencia poblacional que puedes elegir), y lo más interesante de todo por modelo de enfermedad (es capaz de filtrar vcf seleccionando aquellas variantes de novo, heterocigosis recesiva, compuesta, etc.)
  Lo primero que se hace es realizar el filtro por el bed de enriquecimiento que proporciona la casa comercial según el kit de enriquecimiento.
	Según la petición del servicio, se solicita el filtrado de variantes de novo (ONLY include variants at which an offspring has one or two non-inherited alleles AND Exclude variants at which both affected and unaffected subjects have the same heterozygous genotypesc) y en heterocigosis recesiva y compuesta (This function is designed for a disorder suspected to be under compound-heterozygous or recessive inheritance mode, in which both copies of a gene on the two ortholog chromosomes of a patient are damaged by two different mutations or the same mutation. For recessive mode, it simply checks variants with homozygous genotypes in patients. For the compound-heterozygous mode, it can use two different input data, phased genotypes of a patient or unphased genotypes in a trio. Here the trio refers to the two parents and an offspring. When these alleles causing a disease at one locus, it follows the recessive model; and when they are at two loci, it follows the compound-heterozygosity model. In both cases, a gene is hit twice. This is the reason why it has the name 'double-hit gene'.). 
Además de las anotaciones funcionales y predictivas se realizan una serie de filtros: 
Depth < 4
GQ < 10.0
PL < 20
Sequencing quality < 50.0
Population frequency in ANY database (ESP5400,dbsnp141,1kg201305,exac) > 0.005
Los ficheros finales anotados con variantes se pueden encontrar en ..\ANALYSIS\20180220_TRIO04401 \annotation y ..\ANALYSIS\20180220_TRIO04401 \annotation\bedfilter (sin y con filtro por bed de enriquecimiento).
Double-hit variants:..\ANALYSIS\20180220_TRIO04401\annotation\variants_doublehit_excell.xlsx
Novo variants: ..\ANALYSIS\20180220_TRIO04401\annotation\variants_novo_excell.xlsx

La explicación del significado de cada uno de los campos del excell se puede encontrar en: fields_description.xlsx

Estadísticas de las variantes: se ha realizado una serie de estadísticas de las variantes detectadas. Las estadísticas se pueden consultar en el excell: ..\ANALYSIS\20180220_TRIO04401\ stats\vcf_stats.xlxs.
En el enlace se puede leer alguna de las interpretaciones que se pueden hacer con este tipo de estadísticas para el control de calidad del variant calling.
## Annex III

|BAIT_SET|The name of the bait set used in the hybrid selection.|
| --- | --- |
|GENOME_SIZE|The number of bases in the reference genome used for alignment.|
|BAIT_TERRITORY|The number of bases which are localized to one or more baits.|
|TARGET_TERRITORY|The unique number of target bases in the experiment, where the target sequence is usually exons etc.|
|BAIT_DESIGN_EFFICIENCY|The ratio of TARGET_TERRITORY/BAIT_TERRITORY. A value of 1 indicates a perfect design efficiency, while a valud of 0.5 indicates that half of bases within the bait region are not within the target region.|
|TOTAL_READS|The total number of reads in the SAM or BAM file examined.|
|PF_READS|The total number of reads that pass the vendor's filter.|
|PF_UNIQUE_READS|The number of PF reads that are not marked as duplicates.|
|PCT_PF_READS|The fraction of reads passing the vendor's filter, PF_READS/TOTAL_READS.|
|PCT_PF_UQ_READS|The fraction of PF_UNIQUE_READS from the TOTAL_READS, PF_UNIQUE_READS/TOTAL_READS.|
|PF_UQ_READS_ALIGNED|The number of PF_UNIQUE_READS that aligned to the reference genome with a mapping score > 0.|
|PCT_PF_UQ_READS_ALIGNED|The fraction of PF_UQ_READS_ALIGNED from the total number of PF reads.|
|PF_BASES_ALIGNED|The number of PF unique bases that are aligned to the reference genome with mapping scores > 0.|
|PF_UQ_BASES_ALIGNED|The number of bases in the PF_UQ_READS_ALIGNED reads. Accounts for clipping and gaps.|
|ON_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped to the baited regions of the genome.|
|NEAR_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped to within a fixed interval containing a baited region, but not within the baited section per se.|
|OFF_BAIT_BASES|The number of PF_BASES_ALIGNED that are mapped away from any baited region.|
|ON_TARGET_BASES|The number of PF_BASES_ALIGNED that are mapped to a targeted region of the genome.|
|PCT_SELECTED_BASES|The fraction of PF_BASES_ALIGNED located on or near a baited region (ON_BAIT_BASES + NEAR_BAIT_BASES)/PF_BASES_ALIGNED.|
|PCT_OFF_BAIT|The fraction of PF_BASES_ALIGNED that are mapped away from any baited region, OFF_BAIT_BASES/PF_BASES_ALIGNED.|
|ON_BAIT_VS_SELECTED|The fraction of bases on or near baits that are covered by baits, ON_BAIT_BASES/(ON_BAIT_BASES + NEAR_BAIT_BASES).|
|MEAN_BAIT_COVERAGE|The mean coverage of all baits in the experiment.|
|MEAN_TARGET_COVERAGE|The mean coverage of a target region.|
|MEDIAN_TARGET_COVERAGE|The median coverage of a target region.|
|MAX_TARGET_COVERAGE|The maximum coverage of reads that mapped to target regions of an experiment.|
|PCT_USABLE_BASES_ON_BAIT|The number of aligned, de-duped, on-bait bases out of the PF bases available.|
|PCT_USABLE_BASES_ON_TARGET|The number of aligned, de-duped, on-target bases out of all of the PF bases available.|
|FOLD_ENRICHMENT|The fold by which the baited region has been amplified above genomic background.|
|ZERO_CVG_TARGETS_PCT|The fraction of targets that did not reach coverage=1 over any base.|
|PCT_EXC_DUPE|The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.|
|PCT_EXC_MAPQ|The fraction of aligned bases that were filtered out because they were in reads with low mapping quality.|
|PCT_EXC_BASEQ|The fraction of aligned bases that were filtered out because they were of low base quality.|
|PCT_EXC_OVERLAP|The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.|
|PCT_EXC_OFF_TARGET|The fraction of aligned bases that were filtered out because they did not align over a target base.|
|FOLD_80_BASE_PENALTY|The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.|
|PCT_TARGET_BASES_1X|The fraction of all target bases achieving 1X or greater coverage.|
|PCT_TARGET_BASES_2X|The fraction of all target bases achieving 2X or greater coverage.|
|PCT_TARGET_BASES_10X|The fraction of all target bases achieving 10X or greater coverage.|
|PCT_TARGET_BASES_20X|The fraction of all target bases achieving 20X or greater coverage.|
|PCT_TARGET_BASES_30X|The fraction of all target bases achieving 30X or greater coverage.|
|PCT_TARGET_BASES_40X|The fraction of all target bases achieving 40X or greater coverage.|
|PCT_TARGET_BASES_50X|The fraction of all target bases achieving 50X or greater coverage.|
|PCT_TARGET_BASES_100X|The fraction of all target bases achieving 100X or greater coverage.|
|HS_LIBRARY_SIZE|The estimated number of unique molecules in the selected part of the library.|
|HS_PENALTY_10X|The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 10X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.|
|HS_PENALTY_20X|The "hybrid selection penalty" incurred to get 80% of target bases to 20X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.|
|HS_PENALTY_30X|The "hybrid selection penalty" incurred to get 80% of target bases to 30X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.|
|HS_PENALTY_40X|The "hybrid selection penalty" incurred to get 80% of target bases to 40X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 40X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 40 * HS_PENALTY_40X.|
|HS_PENALTY_50X|The "hybrid selection penalty" incurred to get 80% of target bases to 50X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 50X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 50 * HS_PENALTY_50X.|
|HS_PENALTY_100X|The "hybrid selection penalty" incurred to get 80% of target bases to 100X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 100X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 100 * HS_PENALTY_100X.|
|AT_DROPOUT|A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC<=50% regions mapped elsewhere.|
|GC_DROPOUT|A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC>=50% regions mapped elsewhere.|
|HET_SNP_SENSITIVITY|The theoretical HET SNP sensitivity.|
|HET_SNP_Q|The Phred Scaled Q Score of the theoretical HET SNP sensitivity.|

**Bibliografía**
1. *McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., … DePristo, M. a. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. doi:10.1101/gr.107524.110.20*
2. *Li, M.-X., Gui, H.-S., Kwan, J. S. H., Bao, S.-Y., & Sham, P. C. (2012). A comprehensive framework for prioritizing variants in exome sequencing studies of Mendelian diseases. Nucleic acids research, 40(7), e53. doi:10.1093/nar/gkr1257*
