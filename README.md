# **Genomic analysis of Tomato salt stress-induced transcriptome data**

**Author: Li'ang Yu, Andrew Nelson** \
**Date: Aug 31st, 2023**

Tomato salt-induced stress dataset for root tissues contains 60 RNA-seq samples diversity panel includes **60 samples** (**3 biological replicates X 2 conditions X 2 genotypes X 5 time points = 60 samples**) for regulatory studies (See [**meta-info**](https://docs.google.com/spreadsheets/d/1do7XWKwS4BPnlS40jqFlcdd1XHoeCn2S/edit?usp=drive_web&ouid=111454630651995752037&rtpof=true) for sample details). Two species are included along with respective genomes and gene annotations ([**S.lycopersicum**](https://phytozome-next.jgi.doe.gov/info/Slycopersicum_ITAG5_0) and [**S.pimpinellifolium**](https://solgenomics.net/ftp/genomes/Solanum_pimpinellifolium/LA2093/Spimp_LA2093_genome_v1.5/))

The sequencing was performed with Illumina short reads platform (PE 150 bp). Please refer [**Mapping summary**](https://docs.google.com/spreadsheets/d/1OKQqBnini3NFVa3RvF20XS13CJYPJ7Eb/edit#gid=1755113162) for details, including the reads number, Hisat2 mapping rate, FeatureCounts mapping rate, and etc. All raw reads of these accessions are download from [**Novogene**](https://drive.google.com/drive/folders/15ZdFbU9ui28c7eSvuyEwh-sxVHa_TrXn) stored on in-house server workspace (see the following workflow for details). 


- [**Genomic analysis of Tomato salt stress-induced transcriptome data**](#genomic-analysis-of-tomato-salt-stress-induced-transcriptome-data)
  - [Process and quantification of transcriptome data](#process-and-quantification-of-transcriptome-data)
    - [Trim RNA-seq reads using trimmonmatic](#trim-rna-seq-reads-using-trimmonmatic)
    - [**Check the LIB types before alignment**](#check-the-lib-types-before-alignment)
    - [Quantify reads mapping with Hisat2](#quantify-reads-mapping-with-hisat2)
    - [Measure gene expression using FeatureCounts](#measure-gene-expression-using-featurecounts)
      - [S.pimpinellifolium sample quantification](#spimpinellifolium-sample-quantification)
      - [S.lycopersicum sample quantification](#slycopersicum-sample-quantification)
  - [Identify orthologs gene pairs between cultivated and wild species of tomato using liftoff](#identify-orthologs-gene-pairs-between-cultivated-and-wild-species-of-tomato-using-liftoff)
  - [Extract the overlapped genes from two genomes](#extract-the-overlapped-genes-from-two-genomes)
  - [Generate the final report and statistical inforamtion of orthologous genes](#generate-the-final-report-and-statistical-inforamtion-of-orthologous-genes)



## Process and quantification of transcriptome data
### Trim RNA-seq reads using trimmonmatic

```bash
reads_dir="/mnt/Leonof/1_Tomato-SaltStress/01_data"
clean_dir="/mnt/Leonof/1_Tomato-SaltStress/01_data/01_trim"

for i in $(cat $reads_dir/sample.list); 
do 
 trimmomatic PE -threads 20 \
    ${reads_dir}/${i}_1.fq.gz \
    ${reads_dir}/${i}_2.fq.gz \
    ${clean_dir}/${i}_forward_paired.fq.gz ${clean_dir}/${i}_forward_unpaired.fq.gz \
    ${clean_dir}/${i}_reverse_paired.fq.gz ${clean_dir}/${i}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/liangyu/anaconda3/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36
done
```

### **Check the LIB types before alignment**
Input LIBTYPE A to allow Salmon to identify libtype in alignment-based mode. \
The output should be a string like ISR which is stranded data and read comes from the reverse strand.

```bash
fastq_dir="/mnt/Leonof/1_Tomato-SaltStress/01_data/01_trim"
Ref_dir="/mnt/Leonof/1_Tomato-SaltStress/02_Genome"
sample_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data"
output_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon"


###build index for salmon
/home/liangyu/anaconda3/bin/salmon index \
  -t $Ref_dir/Slycopersicum_796_ITAG5.0.fa \
  -i $Ref_dir/Slycopersicum_796_ITAG5.0.fa.salmon

###perform quantification for one sample to check
	/home/liangyu/anaconda3/bin/salmon quant \
    -l A -p 10 \
    -1 $fastq_dir/A1_forward_paired.fq.gz \
    -2 $fastq_dir/A1_reverse_paired.fq.gz \
    -i $Ref_dir/Slycopersicum_796_ITAG5.0.fa.salmon \
    -o /mnt/Leonof/1_Tomato-SaltStress/02_Genome/A1.salmon.count 
```

**See details below**
| TopHat | Paired-end | Single-end |
| ------- | ----------- | ----------- |
| -fr-unstranded |	-l IU |	-l U|
| -fr-firststrand |	-l ISR | -l SR|
| -fr-secondstrand | -l ISF | -l SF|


**Salmon results**
```
{
    "read_files": "[ /mnt/Leonof/1_Tomato-SaltStress/01_data/01_trim/A1_forward_paired.fq.gz, /mnt/Leonof/1_Tomato-SaltStress/01_data/01_trim/A1_reverse_paired.fq.gz]",
    "expected_format": "IU",
    "compatible_fragment_ratio": 1.0,
    "num_compatible_fragments": 14495484,
    "num_assigned_fragments": 14495484,
    "num_frags_with_concordant_consistent_mappings": 12764709,
    "num_frags_with_inconsistent_or_orphan_mappings": 1774854,
    "strand_mapping_bias": 0.5313123863614909,
    "MSF": 0,
    "OSF": 0,
    "ISF": 6782048,
    "MSR": 0,
    "OSR": 0,
    "ISR": 5982661,
    "SF": 901435,
    "SR": 873419,
    "MU": 0,
    "OU": 0,
    "IU": 0,
    "U": 0
}

```
**[2023-09-01 14:12:25.784] [jointLog] [info] Automatically detected most likely library type as IU**

### Quantify reads mapping with Hisat2
```bash
work_dir="/mnt/Leonof/1_Tomato-SaltStress"
clean_dir="/mnt/Leonof/1_Tomato-SaltStress/01_data/01_trim"
bam_dir="/mnt/Leonof/1_Tomato-SaltStress/03_Mapping"

### For domesticated tomato samples 
for i in $(cat $work_dir/Slycopersicum.list); do 

    hisat2 --rna-strandness RF -p 60 --no-softclip \
            --summary-file ${bam_dir}/$i.clean.summary \
            -x /mnt/Leonof/1_Tomato-SaltStress/02_Genome/Slycopersicum_796_ITAG5.0.fa \
            -1 ${clean_dir}/${i}_forward_paired.fq.gz \
    	    -2 ${clean_dir}/${i}_reverse_paired.fq.gz \            
            -S ${bam_dir}/$i.clean.sam

    samtools view -@24 -bS ${bam_dir}/${i}.clean.sam -o ${bam_dir}/${i}.bam
done

### For wild tomato samples 
for i in $(cat $work_dir/Spimpinellifolium.list); do 

    hisat2 --rna-strandness RF -p 60 --no-softclip \
            --summary-file ${bam_dir}/$i.clean.summary \
            -x /mnt/Leonof/1_Tomato-SaltStress/02_Genome/LA2093_genome_v1.5.fa \
            -1 ${clean_dir}/${i}_forward_paired.fq.gz \
 	    -2 ${clean_dir}/${i}_reverse_paired.fq.gz \
            -S ${bam_dir}/$i.clean.sam

    samtools view -@24 -bS ${bam_dir}/${i}.clean.sam -o ${bam_dir}/${i}.bam
done

### clean samples
for i in $(cat ../Spimpinellifolium.list);do 
    mv $i.bam $clean_dir/01_Wild/; 
done

for i in $(cat ../Slycopersicum.list);do
    mv $i.bam $clean_dir/02_Domesticated/; 
done
```

### Measure gene expression using FeatureCounts
#### S.pimpinellifolium sample quantification
```bash
gff_dir="/mnt/Leonof/1_Tomato-SaltStress/02_Genome"
BAM_dir="/mnt/Leonof/1_Tomato-SaltStress/03_Mapping/01_Wild" 
output_dir="/mnt/Leonof/1_Tomato-SaltStress/04_FeatureCounts/01_Wild"

###Perform featureCounts for protein-coding genes (NC-publsihed version)
featureCounts -T 20 -p \
  -a $gff_dir/LA2093_v1.5.gff \
  -o $output_dir/Spimpinellifolium_gene.count \
  -t gene \
  -g ID \
  -s 0 \
  $BAM_dir/*.bam >  S.pimpinellifolium_transcripts_count.log

###Perform featureCounts for protein-coding genes (Fei-lab improved version)
gff_dir="/mnt/Leonof/1_Tomato-SaltStress/02_Genome"
BAM_dir="/mnt/Leonof/1_Tomato-SaltStress/03_Mapping/01_Wild" 
output_dir="/mnt/Leonof/1_Tomato-SaltStress/04_FeatureCounts/01_Wild"

###Perform featureCounts for protein-coding genes
featureCounts -T 20 -p \
  -a $gff_dir/LA2093-rename.gff3 \
  -o $output_dir/Spimpinellifolium_gene.count \
  -t gene \
  -g ID \
  -s 0 \
  $BAM_dir/*.bam >  S.pimpinellifolium_transcripts_count.log

```
**Meta-information summary**
```
### NC published version
||    Features : 35761                                                        ||
||    Meta-features : 35761                                                   ||
||    Chromosomes/contigs : 13  

### Fei lab improved version
||    Features : 36873                                                        ||
||    Meta-features : 36873                                                   ||
||    Chromosomes/contigs : 13     
```

#### S.lycopersicum sample quantification
```bash
gff_dir="/mnt/Leonof/1_Tomato-SaltStress/02_Genome"
BAM_dir="/mnt/Leonof/1_Tomato-SaltStress/03_Mapping/02_Domesticated" 
output_dir="/mnt/Leonof/1_Tomato-SaltStress/04_FeatureCounts/02_Domesticated"

###Perform featureCounts for protein-coding genes
featureCounts -T 20 -p \
  -a $gff_dir/Slycopersicum_796_ITAG5.0.gene_exons.gff3 \
  -o $output_dir/Slycopersicum_gene.count \
  -t gene \
  -g ID \
  -s 0 \
  $BAM_dir/*.bam > S.lycopersicum_transcripts_count.log
```
**Meta-information summary**
```
||    Features : 36648                                                        ||
||    Meta-features : 36648                                                   ||
||    Chromosomes/contigs : 13 
```

## Identify orthologs gene pairs between cultivated and wild species of tomato using liftoff
See [**Liftoff**](https://github.com/agshumate/Liftoff) page for usage details

Liftoff generate fewer genes in soybean samples
```
grep "gene" Gmax_508_Wm82.a4.v1.gene_exons.gff3 | wc -l   52872
grep "gene" Gmax_genes_liftoff.gff3_polished | wc -l   52732

grep "gene" Slyoc_liftoff.gff3_polished | wc -l 38142

```

```bash
lift_dir="/mnt/Leonof/1_Tomato-SaltStress/05_Liftoff"
annotation_dir="/mnt/Leonof/1_Tomato-SaltStress/02_Genome"

# lift S.lycopersicum genes from reference S.lycopersicum to target S.pimpinellifolium
liftoff -polish -copies -sc 0.95 -exclude_partial \
  -s 0.90 -a 0.80 -p 40 \
  -dir /mnt/Leonof/1_Tomato-SaltStress/05_Liftoff/ \
  -o Sslyc_liftoff.gff3 \
  -g /mnt/Leonof/1_Tomato-SaltStress/02_Genome/Slycopersicum_796_ITAG5.0.gene_exons.gff3 \
  /mnt/Leonof/1_Tomato-SaltStress/02_Genome/LA2093_genome_v1.5.fa \
  /mnt/Leonof/1_Tomato-SaltStress/02_Genome/Slycopersicum_796_ITAG5.0.fa

#clean the liftoff gff file

cd $lift_dir/

 awk '$3=="gene" {print}' Sslyc_liftoff.gff3_polished | \
 sed 's/;/\t/g' | awk -F'\t' 'BEGIN {OFS="\t"} {cols=""; for (i=1; i<=NF; i++) {if ($i ~ /coverage/ || $i ~ /sequence/ || $i ~ /valid_ORFs/) {cols = cols (cols=="" ? "" : ",") $i}}; print $9, $1, $4, $5, $7, cols}' | \
 sed 's/,/\t/g' | sed 's/sequence_ID=//g' | sed 's/ID=//g' | sed 's/coverage=//g' \
  > Sslyc_liftoff.gff3_polished.clean

# format liftoff-polished  for bedtools intersect
 awk '{print $2,$3,$4,$1,$5,$6,$7,$8}' Sslyc_liftoff.gff3_polished.clean | sed 's/ /\t/g' >  Sslyc_liftoff.gff3_polished.clean_1

 #format S.pimp gff format (extract the genes only)
 grep 'gene' /mnt/Leonof/1_Tomato-SaltStress/02_Genome/LA2093-rename.gff3 | grep -v "mRNA" | sed 's/;.*//g' | awk '{print $9,$1,$4,$5,$7}' | sed 's/ID=//g' | sed 's/ /\t/g' > Spimp_rename.gff.clean

 # format S.pimp gff  for bedtools intersect
 awk '{print $2,$3,$4,$1,$5}' Spimp_rename.gff.clean | sed 's/ /\t/g' > Spimp_rename.gff.clean_1

 ```
 ## Extract the overlapped genes from two genomes
 ```bash
 
 #identify the overlapped region between lifted region and Spimp gene annoatation
 bedtools intersect -wo -a Spimp_rename.gff.clean_1 \
          -b  Sslyc_liftoff.gff3_polished.clean_1 |\
          sed '1i\Chr Start End S.pimp Strand Chr Start End S.slyc_liftoff Strand Coverage Identity ValidORFs OverlapSize' | \
          awk '{print $4,$1,$2,$3,$5,$9,$6,$7,$8,$10,$11,$12,$13,$14}' | sed 's/ /\t/g' > Spimp_Slyco_genes

# Remove those genes located on different chromosomes
# get S.slyc gene info from S.slyc gff file 
grep -v "#" /mnt/Leonof/1_Tomato-SaltStress/02_Genome/Slycopersicum_796_ITAG5.0.gene_exons.gff3 | sed 's/;/ /g' | \
        awk ' $3=="gene" {print $9,$1,$4,$5,$7}' | sed 's/ /\t/g' | \
sed 's/ID=//g' > Slyco_genes

# Extract the corresponding S.slyc gene infos based on the order of liftoff
for i in $(cat Spimp_Slyco_genes |cut -f6 | cut -f1 -d "_" ); 
do 
  grep $i  Slyco_genes ; 
done | \
sed '1i\Slyco Chr Start End Strand' | sed 's/ /\t/g' >  Slyco_genes.clean

cat Slyco_genes.clean | uniq  > Slyco_genes.final
# combine final outpot 
paste Spimp_Slyco_genes Slyco_genes.clean > Spimp_Slyco_genes_ortholog.txt  

```
## Generate the final report and statistical inforamtion of orthologous genes
```R
setwd('/mnt/Leonof/Xiaodan/soybeans_project/2_analyses/1_liftoff')
ortholog_info <- read.table('Gsoja_Gmax_ortholog.txt', header = T, )
# Count the occurrences of each value
counts <- table(ortholog_info$Gsoja)
counts_column <- counts[ortholog_info$Gsoja]
ortholog <- data.frame(ortholog_info, counts_column)
ortholog$Gsoja_size <- ortholog$End - ortholog$Start
ortholog$Gmax_size <- ortholog$End.2 - ortholog$Start.2
ortholog <-  ortholog[ortholog$Chr == ortholog$Chr.2, ]

pairwise_ortholog <- ortholog[ortholog$Freq == 1, ]
ortholog2 <- subset(ortholog, ortholog$Freq > 1)
# Find the row with the highest value in column2 for each unique value in column1
library(dplyr)
result <- ortholog2 %>%
  arrange(desc(ortholog2$Identity)) %>%
  group_by(ortholog2$Gsoja) %>%
  filter(row_number() == 1)
result <- result[,-24]
ortholog3 <- rbind(pairwise_ortholog, result)

ortholog4 <- ortholog3 %>%
  group_by(ortholog3$Gmax) %>%
  mutate(counts = n())

pairwise_ortholog_1 <- ortholog4[ortholog4$counts == 1, ]
orthologGreat1 <- subset(ortholog4, ortholog4$counts > 1)
result <- orthologGreat1 %>%
  arrange(desc(orthologGreat1$Identity)) %>%
  group_by(orthologGreat1$Gmax) %>%
  filter(row_number() == 1)
result <- result[,-24]

ortholog_final <- rbind(pairwise_ortholog_1, result)
write.table(ortholog_final, file = 'Gsoja_Gmax_ortholog.final.txt', quote = F, col.names = T, row.names = F, sep = '\t')

```


