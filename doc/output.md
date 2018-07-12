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

**Results directory**: ANALYSIS/{ANALYSIS_ID}/01-fastqc
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_fastqc.html`: html report. This file can be opened in your favourite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
   - `{sample_id}/{sample_id}_R[12]_fastqc` : folder with fastqc output in plain text.
   - `{sample_id}/{sample_id}_R[12]_fastqc.zip`: zip compression of above folder.

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used for removal of adapter contamination and trimming of low quality regions. 
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 15 in a 4 nucleotide window.
-  Read lenght < 70

MultiQC reports the percentage of bases removed by trimming in bar plot showing percentage or reads trimmed in forward and reverse.

**Note**:The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the ANALYSIS/{ANALYSIS_ID}/03-preprocQC directory.

**Results directory**: ANALYSIS/{ANALYSIS_ID}/QC/trimmomatic
- There is one folder per sample.
- Files:
   - `{sample_id}/{sample_id}_R[12]_filtered.fastq.gz`: contains high quality reads with both forward and reverse tags surviving.
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
## Picard
Metrics for the analysis of target-capture sequencing experiments are calculated with [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics). The metrics in this class fall broadly into three categories:

- Basic sequencing metrics that are either generated as a baseline against which to evaluate other metrics or because they are used in the calculation of other metrics. This includes things like the genome size, the number of reads, the number of aligned reads etc.
- Metrics that are intended for evaluating the performance of the wet-lab assay that generated the data. This group includes metrics like the number of bases mapping on/off/near baits, %selected, fold 80 base penalty, hs library size and the hs penalty metrics. These metrics are calculated prior to some of the filters are applied (e.g. low mapping quality reads, low base quality bases and bases overlapping in the middle of paired-end reads are all counted).
- Metrics for assessing target coverage as a proxy for how well the data is likely to perform in downstream applications like variant calling. This group includes metrics like mean target coverage, the percentage of bases reaching various coverage levels, and the percentage of bases excluded by various filters. These metrics are computed using the strictest subset of the data, after all filters have been applied.

**Results directory:** ANALYSIS/{ANALYSIS_ID}/99-stats/bamstats
- Files:
   - `hsMetrics_all.out` : summary of the some of the most meaningful columns in picard hsmetrics output for all the samples in the project.
   - `{sample_id}_hsMetrics.out`: full picard hsmetrics output per sample.
   - Description of Picard hsMetrics columns in its output can be found in [AnnexIII](#annex-iii)
   
## Variant Calling
Para el análisis de las variantes germinales en la familia se ha utilizado el workflow de “Best Practices” instaurado por GATK. La documentación completa del framework se puede encontrar (http://www.broadinstitute.org/gatk/guide/best-practices, noviembre 2015). Este workflow se divide en tres pasos diferenciados:

Preprocesamiento de los datos: partiendo del fichero BAM obtenido con BWA y filtrados los duplicados. Se realiza un paso de realineamiento alrededor de los indels para mejorar el alineamiento en zonas donde se sabe que hay problemas (bams realineados en ..\ANALYSIS\20180220_TRIO04401\variants\variants_gatk\realignment) y una recalibración de la calidad de las bases para ajustar las calidades de las bases con un mecanismo de machine learning sabiendo sitios de SNPs conocidos. En la carpeta ..\ANALYSIS\20180220_TRIO04401\variant_calling\variants_gatk\recalibration se encuentran los ficheros pdf (uno por muestra) con gráficas que muestran las estadísticas de calidad antes y después de la recalibración.

Variant Calling: Se realiza la llamada a variantes de forma conjunta para los tres individuos de la familia. El formato saca el genotipo de cada individuo para todas las posiciones en las que alguno de ellos tenga una variante respecto a referencia. Se utiliza el nuevo módulo de variant calling de GATK, HaplotypeCaller que se vale del cálculo de haplotipos para mejorar la llamada de variantes. Además, es capaz de llamar SNPs, indels y algunas variantes estructurales haciendo un ensamblado de novo. 

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

Bibliografía

Li, M.-X., Gui, H.-S., Kwan, J. S. H., Bao, S.-Y., & Sham, P. C. (2012). A comprehensive framework for prioritizing variants in exome sequencing studies of Mendelian diseases. Nucleic acids research, 40(7), e53. doi:10.1093/nar/gkr1257
