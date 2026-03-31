# Gene Regulation in Response to Host Sex and Infection Route in *Brugia pahangi* with New Genome Annotation
Code used *B. pahangi* Differential Expression Analysis and Genome Annotation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19356197.svg)](https://zenodo.org/doi/10.5281/zenodo.19356197)

## Table of Contents
* [*B. pahangi* Genome Annotation](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#b-pahangi-genome-annotation)
    * [Basecall ONT Reads](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#basecall-ont-reads)
        * [Combine FASTQ](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#combine-fastq)
    * [Annotate and Mask Repeates](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#annotate-and-mask-repeats)
    * [Align Reads](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#align-reads)
    * [Remove Duplicates](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#remove-duplicates)
    * [Generate Structural Annotation](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-structural-annotation)
        * [Generate List of Bam Files](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-list-of-bam-files)
        * [Run Braker](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#run-braker)
        * [Combine Annotations](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#combine-annotations)
    * [Generate Updated GFF File](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-updated-gff-file)
    * [Generate Functional Annotation](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-functional-annotation)
    * [Add Locus Tags](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#add-locus-tags)
    * [Check Genbank Compatability](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#check-genbank-compatibility)
* [Pipeline to Generate Counts](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#pipeline-to-generate-counts)
* [Merge and Downsample SQ Samples](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#merge-and-downsample-sq-samples)
    * [Merge SQ Samples](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#merge-sq-samples)
    * [Downsample SQ Samples](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#downsample-sq-samples)
    * [Generate Counts](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-counts)
* [Generate *B. pahangi* GeneInfo](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-b-pahangi-geneinfo)
    * [Generate Polypeptide Files](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-polypeptide-files)
        * [Nuclear Genome](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#nuclear-genome)
        * [Mitochondrial Genome](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#mitochondrial-genome)
    * [Download GO Terms](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#download-go-terms)
    * [Download IPR Terms](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#download-ipr-terms)
    * [Interproscan](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#run-interproscan)
    * [Convert Interproscan to GeneInfo](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#convert-interproscan-to-geneinfo)
* [Generate Updated *B. malayi* GeneInfo](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#generate-updated-b-malayi-geneinfo)
* [Counting *Wolbachia* and Gerbil Reads](https://github.com/christopher-holt/bpahangi_host_sex_effect?tab=readme-ov-file#counting-wolbachia-and-gerbil-reads-per-chromosome)
## *B. pahangi* Genome Annotation
### Basecall ONT Reads
```bash
run_guppy(){

    OUTPUT_DIR=${SAMPLE}
    FAST5_DIR=${OUTPUT_DIR}/fast5
    FASTQ_DIR=${OUTPUT_DIR}/fastq
    
    GUPPY_DIR=${PACKAGE_DIR}/guppy-6.4.2_gpu/
    CONFIG_FILE=${GUPPY_DIR}/data/rna_r9.4.1_70bps_hac.cfg
    GUPPY_BIN_DIR=${GUPPY_DIR}/bin/


    mkdir -p ${FASTQ_DIR}

    ${GUPPY_BIN_DIR}/guppy_basecaller --device cuda:${CUDA} --input_path ${FAST5_DIR} \
    --save_path ${FASTQ_DIR} --config ${CONFIG_FILE} --min_qscore 7 --records_per_fastq 10000000 --gpu_runners_per_device 1 


}

PACKAGE_DIR=
PATH_TO_SAMPLE=
SAMPLE=${PATH_TO_SAMPLE}/20230329-MN21969_BpAF_SQK-RNA002/
CUDA=0
run_guppy

SAMPLE=${PATH_TO_SAMPLE}/20230329-MN23690_BpAM_SQK-RNA002
CUDA=1
run_guppy

```

#### Combine FASTQ
```bash
combine_fastq_files(){
    ONT_FASTQ_DIR=${PATH_TO_SAMPLE}/${SAMPLE}/fastq/pass
    FASTQ=${PATH_TO_SAMPLE}/${SAMPLE}/fastq/${SAMPLE}.fastq
    FASTQ_ZIP=${PATH_TO_SAMPLE}/${SAMPLE}/fastq/${SAMPLE}.fastq.gz

    echo "combining fastq"
    cat ${ONT_FASTQ_DIR}/*fastq > ${FASTQ}
    gzip -c ${FASTQ} > ${FASTQ_ZIP}
    ln -s ${FASTQ_ZIP} ${FASTQ_DIR}/${SAMPLE}.fastq.gz
    echo "fastq combined"
}

SAMPLE=20230329-MN21969_BpAF_SQK-RNA002
combine_fastq_files

SAMPLE=20230329-MN23690_BpAM_SQK-RNA002
combine_fastq_files
```

### Annotate and Mask Repeats
```bash
PACKAGE_DIR=
PROJECT_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/annotate_and_mask_repeats
REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi

REPEATMODELER_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmodeler
REPEATMASKER_SOFTMASKED_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmasker_softmasked
REPEATMASKER_MASKED_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmasker_masked

REPEATMASKER_DIR=${PACKAGE_DIR}/repeatmasker-4.0.7
REPEATMODELER_BIN_DIR=${PACKAGE_DIR}/repeatmodeler-1.0.11/

BPAHANGI_REFERENCE=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta

mkdir -p ${REPEATMODELER_OUTPUT_DIR} ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} ${REPEATMASKER_MASKED_OUTPUT_DIR}

THREADS=10

## Build Database
echo -e "
${REPEATMODELER_BIN_DIR}/BuildDatabase -name ${REPEATMODELER_OUTPUT_DIR}/bpahangi.ncbi.db ${BPAHANGI_REFERENCE} -engine ncbi -dir ${REPEATMODELER_OUTPUT_DIR}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N build_database -l mem_free=30G -wd ${REPEATMODELER_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q

## Model Repeats
echo -e "
${REPEATMODELER_BIN_DIR}/RepeatModeler -pa ${THREADS} -database ${REPEATMODELER_OUTPUT_DIR}/bpahangi.ncbi.db -engine ncbi -dir ${REPEATMODELER_OUTPUT_DIR}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N repeatmodeler -l mem_free=80G -wd ${REPEATMODELER_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid build_database

## Generate Softmasked Fasta
echo -e "
${REPEATMASKER_DIR}/RepeatMasker -pa ${THREADS} -xsmall -gff -e ncbi -lib ${REPEATMODELER_OUTPUT_DIR}/consensi.fa.classified -s -dir ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} ${BPAHANGI_REFERENCE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N soft_repeatmasker -l mem_free=30G -wd ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid repeatmodeler

## Generate Hardmasked Fasta
echo -e "
${REPEATMASKER_DIR}/RepeatMasker -pa ${THREADS} -gff -e ncbi -lib ${REPEATMODELER_OUTPUT_DIR}/consensi.fa.classified -s -dir ${REPEATMASKER_MASKED_OUTPUT_DIR} ${BPAHANGI_REFERENCE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N masked_repeatmasker -l mem_free=30G -wd ${REPEATMASKER_MASKED_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid repeatmodeler
```
## Align Reads
### ONT Reads
```bash
PROJECT_DIR=
PACKAGE_DIR=
MINIMAP2_BIN_DIR=${PACKAGE_DIR}/minimap2-2.17/bin
BAM_DIR=${PROJECT_DIR}/bam
FASTQ_DIR=${PROJECT_DIR}/fastq
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi
ASSEMBLY=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta

mkdir -p ${BAM_DIR}

THREADS=8
map_long_reads_and_calculate_coverage_depth() {
 	echo -e " 
    ${MINIMAP2_BIN_DIR}/minimap2 -t ${THREADS} -ax splice -uf -k14 ${ASSEMBLY} ${FASTQ_DIR}/${SAMPLE}.fastq.gz | \
    ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${BAM_DIR}/${SAMPLE}.sorted.bam && \
        ${SAMTOOLS_BIN_DIR}/samtools index ${BAM_DIR}/${SAMPLE}.sorted.bam &&
       ${SAMTOOLS_BIN_DIR}/samtools flagstat ${BAM_DIR}/${SAMPLE}.sorted.bam > ${BAM_DIR}/${SAMPLE}.flagstat.txt
        " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${BAM_DIR} -N minimap2.${SAMPLE}



}

SAMPLE=20230329-MN21969_BpAF_SQK-RNA002
map_long_reads_and_calculate_coverage_depth

SAMPLE=20230329-MN23690_BpAM_SQK-RNA002
map_long_reads_and_calculate_coverage_depth
```
### Illumina Reads
```bash
PROJECT_DIR=
PACKAGE_DIR=
REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi
HISAT2_DIR=${PACKAGE_DIR}/hisat2-2.1.0
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
REFERENCE_GENOME=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta
INDEXED_REFERENCE=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic
JDK_DIR=${PACKAGE_DIR}/jdk/bin/
TRIMMOMATIC_DIR=${PACKAGE_DIR}/trimmomatic-0.38/

SRR_LIST=${PROJECT_DIR}/annotation.srr.id.list

THREADS=10

FASTQ_DIR=${PROJECT_DIR}/fastq
TRIMMED_READS_DIR=${PROJECT_DIR}/fastq/trimmed_reads
BAM_DIR=${PROJECT_DIR}/bam/

mkdir -p ${BAM_DIR} ${TRIMMED_READS_DIR}

## Create Indexed Reference
${HISAT2_DIR}/hisat2-build ${REFERENCE_GENOME} ${INDEXED_REFERENCE}

for SRR in $(cat ${SRR_LIST} | grep -v "20230329"); do
    SRR_ID=$(echo ${SRR} | cut -d',' -f1)
    SAMPLE_NAME=$(echo ${SRR} | cut -d',' -f3)
    SORTED_BAM_FILE=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam

    FASTQ1=$(echo ${FASTQ_DIR}/${SRR_ID}/${SRR_ID}_R1.fastq.gz)
    FASTQ2=$(echo ${FASTQ_DIR}/${SRR_ID}/${SRR_ID}_R2.fastq.gz)    
    FASTQ1_PAIRED=$(echo ${TRIMMED_READS_DIR}/${SRR_ID}/${SRR_ID}_paired_R1.fastq.gz)
    FASTQ2_PAIRED=$(echo ${TRIMMED_READS_DIR}/${SRR_ID}/${SRR_ID}_paired_R2.fastq.gz)
    FASTQ1_UNPAIRED=$(echo ${TRIMMED_READS_DIR}/${SRR_ID}/${SRR_ID}_unpaired_R1.fastq.gz)
    FASTQ2_UNPAIRED=$(echo ${TRIMMED_READS_DIR}/${SRR_ID}/${SRR_ID}_unpaired_R2.fastq.gz)

    echo -e "
    ${JDK_DIR}/java -Xmx20g -jar ${TRIMMOMATIC_DIR}/trimmomatic-0.38.jar PE -phred33 -trimlog ${TRIMMED_READS_DIR}/trim.log \
    ${FASTQ1} ${FASTQ2} ${FASTQ1_PAIRED} ${FASTQ1_UNPAIRED} ${FASTQ2_PAIRED} ${FASTQ2_UNPAIRED} \
    ILLUMINACLIP:${TRIMMOMATIC_DIR}/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${TRIMMED_READS_DIR} -N trimmomatic.${SRR_ID} -hold_jid fastq.dump.${SRR_ID}


    echo -e "
    ${HISAT2_DIR}/hisat2 -p ${THREADS} --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE} -1 ${FASTQ1_PAIRED} -2 ${FASTQ2_PAIRED} | \
        ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${SORTED_BAM_FILE} && \
       ${SAMTOOLS_BIN_DIR}/samtools index ${SORTED_BAM_FILE}
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread "$THREADS" -q threaded.q -l mem_free=20G -wd "$BAM_DIR" -N hisat2.align.all.${SAMPLE_NAME}

done
```
## Remove Duplicates
```bash
PROJECT_DIR=
PACKAGE_DIR=
THREADS=10

JAVA_BIN_DIR=${PACKAGE_DIR}/jdk/bin
PICARD_DIR=${PACKAGE_DIR}/picard-2.25.3
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin

BAM_DIR=${PROJECT_DIR}/bam/

SRR_LIST=${PROJECT_DIR}/annotation.srr.id.list

for SRR in $(cat ${SRR_LIST} | grep -v "20230329\|Brugia_pahangi_"); do
    SAMPLE_NAME=$(echo ${SRR} | cut -d',' -f3)
    INPUT_BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam
    OUTPUT_BAM=${BAM_DIR}/${SAMPLE_NAME}.dedup.bam
    METRICS_FILE=${BAM_DIR}/${SAMPLE_NAME}.picard.metrics.txt


    echo -e " ${JAVA_BIN_DIR}/java -jar ${PICARD_DIR}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${INPUT_BAM} O=${OUTPUT_BAM} M=${METRICS_FILE} " | \
        qsub -V -q threaded.q -pe thread "$THREADS" -P jhotopp-gcid-proj4b-filariasis -N picard.remove.duplicates.${SAMPLE_NAME} -wd ${BAM_DIR} -l mem_free=70G

    echo -e " ${SAMTOOLS_BIN_DIR}/samtools index -@ ${THREADS} ${OUTPUT_BAM} " | \
        qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N index.picard.${SAMPLE_NAME} -wd ${BAM_DIR} -l mem_free=50G -hold_jid picard.remove.duplicates.${SAMPLE_NAME}


done
```

### Generate Structural Annotation
#### Run Braker
```bash
### Set Directories and File Paths
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

REFERENCE_DIR=${PROJECT_DIR}/annotate_and_mask_repeats/repeatmasker_softmasked
REFERENCE_FASTA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.softmasked.fasta

## Custom Installs
### Genemark
GENEMARK_DIR=${USER_PACKAGE_DIR}/genemark-es-et-4.62 ## Needed to edit perl shebang line for some files


## Set Paths for programs and files
BRAKER_SCRIPT_DIR=${PACKAGE_DIR}/braker-2.1.6/scripts

PROTEIN_OUTPUT_DIR=${PROJECT_DIR}/structural_annotation/braker_protein.final
RNASEQ_OUTPUT_DIR=${PROJECT_DIR}/structural_annotation/braker_rnaseq.final

AUGUSTUS_PKG_DIR=${PACKAGE_DIR}/augustus-3.4.0/
AUGUSTUS_BIN_DIR=${AUGUSTUS_PKG_DIR}/bin/
AUGUSTUS_SCRIPTS_DIR=${AUGUSTUS_PKG_DIR}/scripts/
BAMTOOLS_BIN_DIR=${PACKAGE_DIR}/bamtools-2.5.1/bin/
PYTHON3_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin/
DIAMOND_PATH=${USER_PACKAGE_DIR}/diamond-0.9.24/
CDBFASTA_BIN_DIR=${PACKAGE_DIR}/cdbfasta/
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin/
BLAST_BIN_DIR=${PACKAGE_DIR}/ncbi-blast+-2.8.1/bin/
PROTHINT_BIN_DIR=${GENEMARK_DIR}/ProtHint-2.6.0/bin


SPECIES=Brugiapahangi

export LD_LIBRARY_PATH="${PACKAGE_DIR}/gsl-2.4/lib:$LD_LIBRARY_PATH"

## Comma separated list of the full paths of the RNA-Seq Bam files for braker
BAM_FILES=${PROJECT_DIR}/list_of_bam_files


THREADS=25

## Make sure gm_key is present
if [ ! -f ${HOME}/.gm_key ]; then
 	ln -s ${GENEMARK_DIR}/gm_key ${HOME}/.gm_key
fi

### Set Augustus Dirs with Permissions in Working Dir
mkdir -p ${PROTEIN_OUTPUT_DIR} ${RNASEQ_OUTPUT_DIR}


## RNA-Seq Structrual Annotation
AUGUSTUS_WD_DIR=${RNASEQ_OUTPUT_DIR}/augustus-3.4.0/
AUGUSTUS_CONFIG_PATH=${AUGUSTUS_WD_DIR}/config/
mkdir -p ${AUGUSTUS_WD_DIR}
if [ ! -d ${AUGUSTUS_WD_DIR}/config ]; then

    cp -r ${AUGUSTUS_PKG_DIR}/config ${AUGUSTUS_WD_DIR}
    chmod -R a+rwx ${RNASEQ_OUTPUT_DIR}/augustus-3.4.0

fi

echo -e " 
${BRAKER_SCRIPT_DIR}/braker.pl --species=${SPECIES} \
    --genome=${REFERENCE_FASTA} \
    --bam=$(cat ${BAM_FILES}) \
    --cores=${THREADS} \
    --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
    --AUGUSTUS_BIN_PATH=${AUGUSTUS_BIN_DIR} \
    --AUGUSTUS_SCRIPTS_PATH=${AUGUSTUS_SCRIPTS_DIR} \
    --workingdir=${RNASEQ_OUTPUT_DIR} \
    --BAMTOOLS_PATH=${BAMTOOLS_BIN_DIR} \
    --PYTHON3_PATH=${PYTHON3_BIN_DIR} \
    --GENEMARK_PATH=${GENEMARK_DIR} \
    --DIAMOND_PATH=${DIAMOND_PATH} \
    --CDBTOOLS_PATH=${CDBFASTA_BIN_DIR} \
    --SAMTOOLS_PATH=${SAMTOOLS_BIN_DIR} \
    --useexisting \
    --gff3 \
    --softmasking 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${SPECIES}.rnaseq_braker -l mem_free=100G -wd ${RNASEQ_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q




## Protein Structural Annotation
AUGUSTUS_WD_DIR=${PROTEIN_OUTPUT_DIR}/augustus-3.4.0/
if [ ! -d ${AUGUSTUS_WD_DIR}/config ]; then

    cp -r ${AUGUSTUS_PKG_DIR}/config ${AUGUSTUS_WD_DIR}
    chmod -R a+rwx ${PROTEIN_OUTPUT_DIR}/augustus-3.4.0

fi

echo -e " 
${BRAKER_SCRIPT_DIR}/braker.pl --species=${SPECIES} \
    --genome=${REFERENCE_FASTA} \
    --cores=${THREADS} \
    --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
    --AUGUSTUS_BIN_PATH=${AUGUSTUS_BIN_DIR} \
    --AUGUSTUS_SCRIPTS_PATH=${AUGUSTUS_SCRIPTS_DIR} \
    --workingdir=${PROTEIN_OUTPUT_DIR} \
    --BAMTOOLS_PATH=${BAMTOOLS_BIN_DIR} \
    --PYTHON3_PATH=${PYTHON3_BIN_DIR} \
    --GENEMARK_PATH=${GENEMARK_DIR} \
    --DIAMOND_PATH=${DIAMOND_PATH} \
    --CDBTOOLS_PATH=${CDBFASTA_BIN_DIR}  \
    --BLAST_PATH=${BLAST_BIN_DIR} \
    --SAMTOOLS_PATH=${SAMTOOLS_BIN_DIR} \
    --PROTHINT_PATH=${PROTHINT_BIN_DIR} \
    --prot_seq=${PROJECT_DIR}/databases/proteins.fasta \
    --epmode \
    --useexisting \
    --gff3 \
    --softmasking 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${SPECIES}.protein_braker -l mem_free=50G -wd ${PROTEIN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q

```
#### Combine Annotations
```bash
PROJECT_DIR=
USER_PACKAGE_DIR=

TSEBRA_DIR=${USER_PACKAGE_DIR}/TSEBRA-1.0.3/
ANNOTATION_DIR=${PROJECT_DIR}/structural_annotation
PROT_ONLY=${ANNOTATION_DIR}/braker_protein.final/
RNA_SEQ_ONLY=${ANNOTATION_DIR}/braker_rnaseq.final

${TSEBRA_DIR}/bin/tsebra.py -g ${RNA_SEQ_ONLY}/augustus.hints.gtf,${PROT_ONLY}/augustus.hints.gtf -c ${TSEBRA_DIR}/config/default.cfg -e ${RNA_SEQ_ONLY}/hintsfile.gff,${PROT_ONLY}/hintsfile.gff -o ${ANNOTATION_DIR}/aviteae_structural_annotation.final.gtf

```

### Generate Updated GFF File
```bash
SCRIPT_DIR=
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/structural_annotation/
BPAHANGI_AA=${ANNOTATION_DIR}/bpahangi_structural_annotation.final.aa
GTF=${ANNOTATION_DIR}/bpahangi_structural_annotation.final.gtf
GFF=${ANNOTATION_DIR}/bpahangi_structural_annotation.final.gff3


REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi
BPAHANGI_REFERENCE=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta

THREADS=1

echo -e "
${PACKAGE_DIR}/r-4.1.2/bin/Rscript ${SCRIPT_DIR}/convert_tsebra_gtf_to_gff.R --gtf ${GTF} --gff ${GFF}
" | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N convert_file_format -wd ${ANNOTATION_DIR} -l mem_free=50G



echo -e "
${USER_PACKAGE_DIR}/gffread-0.12.7/gffread -y ${BPAHANGI_AA} -g ${BPAHANGI_REFERENCE} ${GFF}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${ANNOTATION_DIR} -N gff_read -hold_jid convert_file_format

```

### Generate Functional Annotation
```bash
SCRIPT_DIR=
PROJECT_DIR=

AA_HINTS=${PROJECT_DIR}/structural_annotation/bpahangi_structural_annotation.final.aa
BPAHANGI_ANNOTATION=${PROJECT_DIR}/structural_annotation/bpahangi_structural_annotation.final.gff3
FUNCTIONAL_ANNOTATION_DIR=${PROJECT_DIR}/functional_annotation/tsebra.final
BPAHANGI_FASTA=${PROJECT_DIR}/reference/brugia_pahangi/GCA_012070555.1_ASM1207055v1_genomic.fasta


sh ${SCRIPT_DIR}/tsebra_functional_annotation.sh \
    -a ${AA_HINTS} \
    -s ${BPAHANGI_ANNOTATION} \
    -o ${FUNCTIONAL_ANNOTATION_DIR} \
    -f ${BPAHANGI_FASTA} \
    -g ${FUNCTIONAL_ANNOTATION_DIR}/bpahangi_functional_annotation.final.gff3

```

### Add Locus Tags
The QQG55 tag was provided by NCBI
```bash
PROJECT_DIR=
FUNCTIONAL_ANNOTATION_DIR=${PROJECT_DIR}/functional_annotation/tsebra.final
## BIOCODE_DIR obtained from https://github.com/jorvis/biocode
BIOCODE_DIR=

${BIOCODE_DIR}/gff/add_gff3_locus_tags.py -i ${FUNCTIONAL_ANNOTATION_DIR}/bpahangi_functional_annotation.final.gff3 -o ${FUNCTIONAL_ANNOTATION_DIR}/bpahangi_functional_annotation.final.with_locus_tags.gff3 -p QQG55 -a 4 -n 5 -s 5
```
### Check Genbank Compatibility
```bash
PROJECT_DIR=
GENBANK_DIR=${PROJECT_DIR}/genbank_submission

## Nuclear Genome
GFF3=${GENBANK_DIR}/bpahangi_functional_annotation.final.with_locus_tags.gff3 ## Symlink of ${FUNCTIONAL_ANNOTATION_DIR}/bpahangi_functional_annotation.final.with_locus_tags.gff3
SBT=${GENBANK_DIR}/Bpahangi.sbt
CMT=${GENBANK_DIR}/Bpahangi.cmt
FSA=${GENBANK_DIR}/Bpahangi_nucgenome_JAAVKF.fasta ## Genome assembly file excluding the mitochondria

${USER_PACKAGE_DIR}/tbl2asn/linux64.table2asn_GFF -J -c w -t ${SBT} -w ${CMT} -i ${FSA} -f ${GFF3} -j "[gcode=1][organism=Brugia pahangi][strain=FR3][host=Meriones unguiculatus][country=USA][collection-date=24-Jul-2017]" -M n -V vb -Z -euk


## Mitochondria
GFF3=${GENBANK_DIR}/bpahangi_functional_annotation.final.with_locus_tags.mitochondria.gff3 ## Manual annotation of B. pahangi mitochondria
FSA=${GENBANK_DIR}/Bpahangi_mitogenome_JAAVKF.fasta ## Genome assembly file including only the mitochondria

${USER_PACKAGE_DIR}/tbl2asn/linux64.table2asn_GFF -J -c w -t ${SBT} -w ${CMT} -i ${FSA} -f ${GFF3} -j "[mgcode=5][organism=Brugia pahangi][strain=FR3][host=Meriones unguiculatus][country=USA][collection-date=24-Jul-2017]" -M n -V vb -Z -euk

```



## Pipeline to Generate Counts
```bash
PACKAGE_DIR=
USER_PACKAGE_DIR=
PROJECT_DIR=

SRR_LIST=${PROJECT_DIR}/srr.id.list
REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi

## Reference Genome and GFF3
REFERENCE_GENOME=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta
GFF3=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.gff
INDEXED_REFERENCE=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic

## Set up variables (programs and directories)
SRA_TOOLKIT_BIN_DIR=${PACKAGE_DIR}/sratoolkit-2.10.9/bin/
JDK_DIR=${PACKAGE_DIR}/jdk/bin/
TRIMMOMATIC_DIR=${PACKAGE_DIR}/trimmomatic-0.38/
HISAT2_DIR=${PACKAGE_DIR}/hisat2-2.2.1/
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
PICARD_DIR=${PACKAGE_DIR}/picard-2.25.3
PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin

THREADS=6

## Set Directories and Output Files
SRR_DIR=${PROJECT_DIR}/fastq/original_reads
BPAHANGI_TRIMMED=${PROJECT_DIR}/fastq/trimmed_reads
BAM_DIR=${PROJECT_DIR}/bam
COUNTS_DIR=${PROJECT_DIR}/counts
COMBINED_COUNTS=${COUNTS_DIR}/combined.final.counts

mkdir -p ${SRR_DIR} ${BPAHANGI_TRIMMED} ${BAM_DIR} ${COUNTS_DIR}

## Create Index of Reference Genome
echo -e "
${HISAT2_DIR}/hisat2-build ${REFERENCE_GENOME} ${INDEXED_REFERENCE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=20G -wd ${REFERENCE_DIR} -N hisat2.build



for SRR in $(cat ${SRR_LIST}); do
    SRR_ID=$(echo ${SRR} | cut -d',' -f1)
    SAMPLE=$(echo ${SRR} | cut -d',' -f3)
    SORTED_BAM_FILE=${BAM_DIR}/${SAMPLE}.sorted.bam


    ## Trimmed Reads Files
    TRIMMED_READS_DIR=${BPAHANGI_TRIMMED}/${SRR_ID}
    FASTQ1=${SRR_DIR}/${SRR_ID}/${SRR_ID}_1.fastq.gz
    FASTQ2=${SRR_DIR}/${SRR_ID}/${SRR_ID}_2.fastq.gz
    FASTQ1_PAIRED=${TRIMMED_READS_DIR}/${SRR_ID}_paired_1.fastq.gz
    FASTQ2_PAIRED=${TRIMMED_READS_DIR}/${SRR_ID}_paired_2.fastq.gz
    FASTQ1_UNPAIRED=${TRIMMED_READS_DIR}/${SRR_ID}_unpaired_1.fastq.gz
    FASTQ2_UNPAIRED=${TRIMMED_READS_DIR}/${SRR_ID}_unpaired_2.fastq.gz

    mkdir -p ${TRIMMED_READS_DIR}

    ## Dedup Bam Files
    INPUT_BAM=${BAM_DIR}/${SAMPLE}.sorted.bam
    DEDUP_BAM=${BAM_DIR}/${SAMPLE}.dedup.bam
    METRICS_FILE=${BAM_DIR}/${SAMPLE}.picard.metrics.txt

    ## 1. Download .sra file using pre-fetch
    echo -e "
    ${SRA_TOOLKIT_BIN_DIR}/prefetch ${SRR_ID} --max-size 100G --output-directory ${SRR_DIR}
    " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=20G -wd ${SRR_DIR} -N prefetch.${SRR_ID}

    ## 2. Unpack .sra file into fastq.gz files
    echo -e "
    ${SRA_TOOLKIT_BIN_DIR}/fastq-dump --split-files ${SRR_ID} --gzip -O ${SRR_DIR}/${SRR_ID}
    " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=20G -wd ${SRR_DIR} -N fastq.dump.${SRR_ID} -hold_jid prefetch.${SRR_ID}

    ## 3. Trim Reads
    echo -e "
    ${JDK_DIR}/java -Xmx20g -jar ${TRIMMOMATIC_DIR}/trimmomatic-0.38.jar PE -phred33 -trimlog ${TRIMMED_READS_DIR}/trim.log \
    ${FASTQ1} ${FASTQ2} ${FASTQ1_PAIRED} ${FASTQ1_UNPAIRED} ${FASTQ2_PAIRED}  \
    ${FASTQ2_UNPAIRED} \
    ILLUMINACLIP:${TRIMMOMATIC_DIR}/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${TRIMMED_READS_DIR} -N trimmomatic.${SRR_ID} -hold_jid fastq.dump.${SRR_ID}

    ## 4. Align Reads to Reference
    echo -e "
    ${HISAT2_DIR}/hisat2 --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE} -1 ${FASTQ1_PAIRED} -2 ${FASTQ2_PAIRED} | \
    ${SAMTOOLS_BIN_DIR}/samtools sort -o ${SORTED_BAM_FILE} && \
    ${SAMTOOLS_BIN_DIR}/samtools index ${SORTED_BAM_FILE}
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=70G -wd ${BAM_DIR} -N hisat2.align.${SAMPLE} -hold_jid trimmomatic.${SRR_ID},hisat2.build


    ## 5. Remove Duplicates
    echo -e "
    ${JDK_DIR}/java -Xmx20g -jar ${PICARD_DIR}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${INPUT_BAM} O=${DEDUP_BAM} M=${METRICS_FILE}
    " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N picard.remove.duplicates.${SAMPLE} -wd ${BAM_DIR} -l mem_free=70G -hold_jid hisat2.align.${SAMPLE}

    echo -e "
    ${SAMTOOLS_BIN_DIR}/samtools index ${DEDUP_BAM}
    " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N index.picard.${SAMPLE} -wd ${BAM_DIR} -l mem_free=50G -hold_jid picard.remove.duplicates.${SAMPLE}


    ## 6. Generate Counts
    echo -e "
    ${PYTHON_BIN_DIR}/htseq-count -n ${THREADS} -s reverse --max-reads-in-buffer 3000000000 -r pos --nonunique none -f bam -m union -t gene --idattr ID ${DEDUP_BAM} ${GFF3} \
    | awk -v a=${SAMPLE} '{print \$0, a}' | sed -e 's/ /\t/g' >> ${COMBINED_COUNTS}
    " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N htseq_counts.${SAMPLE} -wd ${COUNTS_DIR} -l mem_free=50G -hold_jid index.picard.${SAMPLE}

done

```
## Merge and Downsample SQ Samples
### Merge SQ Samples
```bash
PACKAGE_DIR=
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.20/bin
cd ${BAM_DIR}

sort_and_index_bam_file(){
    ${SAMTOOLS_BIN_DIR}$/samtools sort -@ 8 -o ${SORTED_MERGED_DEDUP_BAM} ${MERGED_DEDUP_BAM}
    ${SAMTOOLS_BIN_DIR}$/samtools index ${SORTED_MERGED_DEDUP_BAM}
}

## Male worms from SQ Male Gerbils
MERGED_DEDUP_BAM=M_SQ_Bp-males.dedup.bam
SORTED_MERGED_DEDUP_BAM=M_SQ_Bp-males.sorted.dedup.bam
${SAMTOOLS_BIN_DIR}$/samtools merge -o ${MERGED_DEDUP_BAM} -f -@ 8 M055_SQ_Bp-males.dedup.bam M057_SQ_Bp-males.dedup.bam M069_SQ_Bp-males.dedup.bam
sort_and_index_bam_file

## Female worms from SQ Male Gerbils
MERGED_DEDUP_BAM=M_SQ_Bp-females.dedup.bam
SORTED_MERGED_DEDUP_BAM=M_SQ_Bp-females.sorted.dedup.bam
${SAMTOOLS_BIN_DIR}$/samtools merge -o ${MERGED_DEDUP_BAM}  -f -@ 8 M055_SQ_Bp-females.dedup.bam M057_SQ_Bp-females.dedup.bam M069_SQ_Bp-females.dedup.bam
sort_and_index_bam_file

## Male worms from SQ Female Gerbils
MERGED_DEDUP_BAM=F_SQ_Bp-males.dedup.bam
SORTED_MERGED_DEDUP_BAM=F_SQ_Bp-males.sorted.dedup.bam
${SAMTOOLS_BIN_DIR}$/samtools merge -o ${MERGED_DEDUP_BAM} -f -@ 8 F024_SQ_Bp-males.dedup.bam F022_SQ_Bp-males.dedup.bam F015_SQ_Bp-males.dedup.bam
sort_and_index_bam_file

## Female worms from SQ Female Gerbils
MERGED_DEDUP_BAM=F_SQ_Bp-females.dedup.bam
SORTED_MERGED_DEDUP_BAM=F_SQ_Bp-females.sorted.dedup.bam
sort_and_index_bam_file
```
### Downsample SQ Samples
```bash
PROJECT_DIR=
PACKAGE_DIR=
JDK_DIR=${PACKAGE_DIR}/jdk/bin/
PICARD_DIR=${PACKAGE_DIR}/picard-2.25.3
BAM_DIR=${PROJECT_DIR}/bam

cd ${BAM_DIR}

downsample_bam(){
    ${JDK_DIR}/java -jar ${PICARD_DIR}/picard.jar DownsampleSam I=${SORTED_MERGED_DEDUP_BAM} O=${DOWNSAMPLED_BAM} P=0.6
}

## Male worms from SQ Male Gerbils
SORTED_MERGED_DEDUP_BAM=M_SQ_Bp-males.sorted.dedup.bam
DOWNSAMPLED_BAM=M_SQ_Bp-males.sorted.dedup.100m.bam
downsample_bam

## Female worms from SQ Male Gerbils
SORTED_MERGED_DEDUP_BAM=M_SQ_Bp-females.sorted.dedup.bam
DOWNSAMPLED_BAM=M_SQ_Bp-females.sorted.dedup.100m.bam
downsample_bam

## Female worms from SQ Female Gerbils
SORTED_MERGED_DEDUP_BAM=F_SQ_Bp-females.sorted.dedup.bam
DOWNSAMPLED_BAM=F_SQ_Bp-females.sorted.dedup.100m.bam
downsample_bam

## Male worms from SQ Female Gerbils
SORTED_MERGED_DEDUP_BAM=F_SQ_Bp-males.sorted.dedup.bam
DOWNSAMPLED_BAM=F_SQ_Bp-males.sorted.dedup.100m.bam
downsample_bam
```
### Generate Counts
```bash
PROJECT_DIR=
PACKAGE_DIR=
REFERENCE_DIR=${PROJECT_DIR}/reference/brugia_pahangi
GFF3=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.gff
PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin
COUNTS_DIR=${PROJECT_DIR}/counts
THREADS=16

generate_counts(){
    COUNTS_FILE=${COUNTS_DIR}/${SAMPLE_NAME}.counts
    DOWNSAMPLED_BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.dedup.100m.bam
    ${PYTHON_BIN_DIR}/htseq-count -n ${THREADS} -s reverse --max-reads-in-buffer 3000000000 -r pos --nonunique none -f bam -m union -t gene --idattr ID ${DOWNSAMPLED_BAM} ${GFF3} | awk -v a=${SAMPLE_NAME} '{print $0, a}' | sed -e 's/ /\t/g' > ${COUNTS_FILE}
}

## Male worms from SQ Male Gerbils
SAMPLE_NAME=M_SQ_Bp-males
generate_counts

## Female worms from SQ Male Gerbils
SAMPLE_NAME=M_SQ_Bp-females
generate_counts

## Female worms from SQ Female Gerbils
SAMPLE_NAME=F_SQ_Bp-females
generate_counts

## Male worms from SQ Female Gerbils
SAMPLE_NAME=F_SQ_Bp-males
generate_counts
```


## Generate *B. pahangi* GeneInfo
### Generate Polypeptide Files
#### Nuclear Genome
```bash
## 1. Create fasta of the protein sequences for the nuclear genome

THREADS=4

NUCLEAR_AA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_no_mitochondria.aa
NUCLEAR_FNA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_no_mitochondria.fna
NUCLEAR_GFF=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_no_mitochondria.gff
REFERENCE_GENOME=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.fasta


GFFREAD_DIR=${USER_PACKAGE_DIR}/gffread-0.12.7
EMBOSS_BIN_DIR=${PACKAGE_DIR}/emboss-6.6.0/bin

## Remove Mitochondria from GFF3
cat ${GFF3} | grep -v "CM022469.1" > ${NUCLEAR_GFF3}

## Generate FNA file
${GFFREAD_DIR}/gffread -x ${NUCLEAR_FNA} -g ${REFERENCE_GENOME} ${NUCLEAR_GFF}
## Convert FNA to Amino Acid
${EMBOSS_BIN_DIR}/transeq -sequence ${NUCLEAR_FNA} -outseq ${NUCLEAR_AA} -table 1 -clean
```
#### Mitochondrial Genome
```bash
## 2. Create fasta of the protein sequences for the mitochondrial genome

MITO_AA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_mitochondria.aa
MITO_FNA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_mitochondria.fna
MITO_GFF=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic_mitochondria.gff

GFFREAD_DIR=${USER_PACKAGE_DIR}/gffread-0.12.7
EMBOSS_BIN_DIR=${PACKAGE_DIR}/emboss-6.6.0/bin

## Select only Mitochondria
cat ${GFF3} | grep "CM022469.1" > ${MITO_GFF}

## Generate FNA file
${GFFREAD_DIR}/gffread -x ${MITO_FNA} -g ${REFERENCE_GENOME} ${MITO_GFF}
## Convert FNA to Amino Acid
${EMBOSS_BIN_DIR}/transeq -sequence ${MITO_FNA} -outseq ${MITO_AA} -table 5 -clean
```
### Download GO Terms
```bash
## 3. Download GO Terms
R_BIN_DIR=${PACKAGE_DIR}/r-4.1.2/bin
${R_BIN_DIR}/Rscript ${SCRIPTS}/download_go_terms.R --output ${PROJECT_DIR}/goid_and_description.map
```
### Download IPR Terms
```bash
## 4. Download IPR Terms
PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin

## Code is available at https://www.ebi.ac.uk/interpro/result/download/#/entry/InterPro/|tsv

${PYTHON_BIN_DIR}/python download_ipr_terms.py > ${PROJECT_DIR}/interproid_and_description.map
```

### Run InterproScan
```bash
# 5. Run Interproscan
INTERPROSCAN_DIR=${USER_PACKAGE_DIR}/interproscan-5.56-89.0
OUTPUT_DIR=${PROJECT_DIR}/interproscan
BPAHANGI_AA=${REFERENCE_DIR}/GCA_012070555.1_ASM1207055v1_genomic.aa

cat ${NUCLEAR_AA} ${MITO_AA} > ${BPAHANGI_AA}


mkdir -p ${OUTPUT_DIR}

export PATH=${PACKAGE_DIR}/jdk-18/bin:$PATH


if [ ! -e ${HOME}/lib/libpcre.so ]; then
    ln -s /usr/lib64/libpcre.so.1.2.10 ${HOME}/lib/libpcre.so
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${HOME}/lib

THREADS=15

echo -e "
${INTERPROSCAN_DIR}/interproscan.sh -cpu ${THREADS} -i ${BPAHANGI_AA}  -d ${OUTPUT_DIR} -f tsv,gff3 -goterms -iprlookup -verbose
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=40G -wd ${OUTPUT_DIR} -N interpro

```
### Convert InterproScan to GeneInfo
```bash
## 6. Convert Interproscan to geneinfo
PROJECT_DIR=
OUTPUT_DIR=${PROJECT_DIR}/interproscan
INTERPROSCAN_TSV=${OUTPUT_DIR}/GCA_012070555.1_ASM1207055v1_genomic.aa.tsv
INTERPROSCAN_FORMATTED=${OUTPUT_DIR}/GCA_012070555.1_ASM1207055v1_genomic.aa.formatted.tsv
GOMAP=${PROJECT_DIR}/goid_and_description.map
IPRMAP=${PROJECT_DIR}/interproid_and_description.map
COUNTS_DIR=${PROJECT_DIR}/counts
COMBINED_COUNTS=${COUNTS_DIR}/combined.final.counts

### Remove _1 added by interproscan, remove -2, swap rna- to gene-
sed 's/_1*//2' ${INTERPROSCAN_TSV} | sed -e 's/-[0-9]//1' | sed -e 's/rna-QQG/gene-QQG/1' > ${INTERPROSCAN_FORMATTED}


THREADS=1
## Convert interproscan to gene.info file
echo -e "
${SCRIPTS}/interproscan_to_geneinfo.R --interproscan ${INTERPROSCAN_FORMATTED}  --counts ${COMBINED_COUNTS} --gomap ${GOMAP} --iprmap ${IPRMAP} --out ${OUTPUT_DIR}/bpahangi_gene.info
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${OUTPUT_DIR} -N interproscan_to_geneinfo

```

## Generate Updated *B. malayi* GeneInfo
```bash
PROJECT_DIR=
SCRATCH_DIR=
INTERPROSCAN_DIR=${USER_PACKAGE_DIR}/interproscan-5.56-89.0
OUTPUT_DIR=${PROJECT_DIR}/bmalayi_interproscan

BMALAYI_AA=${REFERENCE_DIR}/b_malayi.PRJNA10729.WS276.genomic.aa ## Created using GFFRead and Emboss Transeq, same as above

mkdir -p ${OUTPUT_DIR}

export PATH=${PACKAGE_DIR}/jdk-18/bin:$PATH


if [ ! -e ${HOME}/lib/libpcre.so ]; then
    ln -s /usr/lib64/libpcre.so.1.2.10 ${HOME}/lib/libpcre.so
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${HOME}/lib

${INTERPROSCAN_DIR}/interproscan.sh \
    -cpu ${THREADS} \
    -i ${BMALAYI_AA}  \
    -d ${OUTPUT_DIR} \
    -f tsv,gff3 -goterms -iprlookup -verbose \
    --tempdir ${SCRATCH_DIR} \
    --applications Pfam,SMART,CDD,TIGRFAM,SuperFamily


sed -e 's/_1//1' ${OUTPUT_DIR}/b_malayi.PRJNA10729.WS276.genomic.aa.tsv | sed -e 's/-[0-9]//1' > ${OUTPUT_DIR}/b_malayi.PRJNA10729.WS276.genomic.aa.formatted.tsv

GOMAP=${HOME}/maps/goid_and_description.map
IPRMAP=${HOME}/maps/interproid_and_description.map
GENE_NAMES=${OUTPUT_DIR}/WB_terms.txt ## List of all WB gene names

## b_malayi.PRJNA10729.WS276.genomic.aa.formatted.names_updated.tsv was created by swapping out the Bm names for the corresponding WB terms. 
## The Bm_to_WB.tsv table is available under input files and was created by scripts/convert_Bm_to_WB_interproscan.R

${SCRIPTS}/interproscan_to_geneinfo.R \
    --interproscan ${OUTPUT_DIR}/b_malayi.PRJNA10729.WS276.genomic.aa.formatted.names_updated.tsv  \
    --counts ${GENE_NAMES} \
    --gomap ${GOMAP} \
    --iprmap ${IPRMAP} \
    --out ${OUTPUT_DIR}/bmalayi_gene.info
```

## Counting *Wolbachia* and Gerbil Reads Per Chromosome
```bash
## Generated by hisat2-build, used a composite reference of Brugia pahangi + mitochondria (GCA_012070555) and it's Wolbachia (GCA_012030695)
INDEXED_REFERENCE=${REFERENCE_DIR}/GCA_bpahangi_wBp_combined
SPECIES=brugia_pahangi
## To get statistics for the gerbil genome, this reference was used
INDEXED_REFERENCE=${REFERENCE_DIR}/GCF_030254825.1_Bangor_MerUng_6.1_genomic
SPECIES=gerbil
```

```bash
PROJECT_DIR=
PACKAGE_DIR=
HISAT2_DIR=${PACKAGE_DIR}/hisat2-2.2.1/
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
PICARD_DIR=${PACKAGE_DIR}/picard-2.25.3
JDK_DIR=${PACKAGE_DIR}/jdk/bin/

SRR_LIST=${PROJECT_DIR}/groups.txt

SAMPLES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SRR_LIST})
SRR_ID=$(echo ${SAMPLES} | cut -d',' -f1)
SAMPLE_NAME=$(echo ${SAMPLES} | cut -d',' -f3)

## Trimmed reads
FASTQ1_PAIRED=${PROJECT_DIR}/trimmed_reads/${SRR_ID}_paired_1.fastq.gz
FASTQ2_PAIRED=${PROJECT_DIR}/trimmed_reads/${SRR_ID}_paired_2.fastq.gz

## Output bam files, ${SPECIES} defined above as either gerbil or brugia pahangi
BAM_DIR=${PROJECT_DIR}/bam_${SPECIES}
SORTED_BAM_FILE=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam
DEDUP_BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.dedup.bam
METRICS_FILE=${BAM_DIR}/${SAMPLE_NAME}.picard.metrics.txt

## Set idxstats output files
IDX_DIR=${PROJECT_DIR}/idxstats_${SPECIES}
IDXSTAT_OUT=${IDX_DIR}/${SAMPLE_NAME}.idxstat.txt

mkdir ${BAM_DIR} ${IDX_DIR}

## Map trimmed reads to the genome, indexed genome defined above
${HISAT2_DIR}/hisat2 --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE} -1 ${FASTQ1_PAIRED} -2 ${FASTQ2_PAIRED} | \
    ${SAMTOOLS_BIN_DIR}/samtools sort -o ${SORTED_BAM_FILE}

## Index Bam File
${SAMTOOLS_BIN_DIR}/samtools index ${SORTED_BAM_FILE}

## Deduplicate Bam Files
${JDK_DIR}/java -Xmx20g -jar ${PICARD_DIR}/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -I ${SORTED_BAM_FILE} -O ${DEDUP_BAM} -M ${METRICS_FILE}

## Index Deduplicated Bam Files
${SAMTOOLS_BIN_DIR}/samtools index ${DEDUP_BAM}

## Count the number of reads per chromosome
${SAMTOOLS_BIN_DIR}/samtools idxstats -@ ${THREADS} ${DEDUP_BAM} > ${IDXSTAT_OUT}
```
