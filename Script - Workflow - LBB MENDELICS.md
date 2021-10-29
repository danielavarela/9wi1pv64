##Script - Workflow - LBB MENDELICS

#####PREPARO DO AMBIENTE

```bash
apt-get uptade install tree fastqc bwa samtools bedtools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
unzip gatk.4.2.2.0
rm gatk.4.2.2.0.zip
wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar
find / -name bwa
```

##### DOWNLOAD DOS DADOS EM FORMATO FASTQ

** Baixei as amostra-lbb_R1.fq.gz amostra-lbb_R2.fq.gz e subi para o JupyterLab **

```bash
zcat amostra-lbb_R1.fq.gz | head -n 4
```
##### PROCESSAMENTO DOS DADOS

```bash
mkdir dados
mkdir dados/fastq
mkdir dados/bwa
mkdir dados/picard
mkdir dados/fastqc
mkdir dados/bedtools
mkdir dados/annovar
mkdir dados/gatk

mv amostra-lbb_R1.fq.gz ./dados/fastq
mv amostra-lbb_R2.fq.gz ./dados/fastq

tree dados
```
##### CONTROLE DE QUALIDADE DOS READS

```bash
time fastqc -o dados/fastqc \
  dados/fastq/amostra-lbb_R1.fq.gz 
time fastqc -o  dados/fastqc \
    dados/fastq/amostra-lbb_R2.fq.gz 

tree dados
```
_Podemos usar o programa cutadapt para remover reads muito curtos, muito longos ou com adaptadores_

```jupyter
pip install --upgrade cutadapt
```
```bash
time cutadapt --minimum-length 90 --maximum-length 151 \
  -o dados/fastq/amostra-lbb_R1_cutadapt.fastq \
  dados/fastq/amostra-lbb_R1.fq.gz 
time cutadapt --minimum-length 90 --maximum-length 151 \
  -o dados/fastq/amostra-lbb_R2_cutadapt.fastq \
  dados/fastq/amostra-lbb_R2.fq.gz 

tree dados
```
_Vamos gerar os relatórios de qualidade novamente, agora com os reads trimados_

```bash
time fastqc -o dados/fastqc \
  dados/fastq/amostra-lbb_R1_cutadapt.fastq
time fastqc -o dados/fastqc \
  dados/fastq/amostra-lbb_R2_cutadapt.fastq
```
##### BAIXANDO E INDEXANDO O GENOMA DE REFERÊNCIA

```bash
mkdir referencia
mkdir referencia/grch38

#Baixei o arquivo fasta disponivel de sequencia referência com o cromossomo 22 do genoma humano GRCh38 e subi para o JupyterLab
 (Descompactar o arquivo gz do cromossomo 22)
   gunzip grch38.chr22.fasta

time bwa index -a bwtsw referencia/grch38/grch38.chr22.fasta

tree dados referencia

#criar o indice da referencia (chr22)
time samtools faidx referencia/grch38/grch38.chr22.fasta

java -jar picard.jar CreateSequenceDictionary \
REFERENCE=referencia/grch38/grch38.chr22.fasta \
OUTPUT=referencia/grch38/grch38.chr22.dict

tree dados referencia
```
##### ALINHAMENTO NO GENOMA DE REFERÊNCIA

```bash
bwa

#Etapa de mapeamento dos reads em relação a referencia do chr22
NOME=DANI
Biblioteca=FocusDNA
Plataforma=Illumina

time bwa mem -M -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" \
referencia/grch38/grch38.chr22.fasta \
dados/fastq/amostra-lbb_R1.fq.gz \
dados/fastq/amostra-lbb_R2.fq.gz dados/bwa/AMOSTRA01.sam 

head dados/bwa/AMOSTRA01.sam

time samtools fixmate dados/bwa/AMOSTRA01.sam dados/bwa/AMOSTRA01.bam
time samtools sort -O bam -o dados/bwa/AMOSTRA01_sorted.bam dados/bwa/AMOSTRA01.bam
time samtools index dados/bwa/AMOSTRA01_sorted.bam

time samtools view dados/bwa/AMOSTRA01_sorted.bam | head

tree dados
```
##### CONFERIR COBERTURA

```bash
#ls -lh dados/bwa/AMOSTRA01_sorted.bam
bedtools bamtobed -i dados/bwa/AMOSTRA01_sorted.bam \
>dados/bedtools/AMOSTRA01_sorted.bed

head dados/bedtools/AMOSTRA01_sorted.bed

bedtools merge -i dados/bedtools/AMOSTRA01_sorted.bed >dados/bedtools/AMOSTRA01_merged.bed
bedtools sort -i dados/bedtools/AMOSTRA01_merged.bed >dados/bedtools/AMOSTRA01_merged_sorted.bed

wc -l dados/bedtools/AMOSTRA01_sorted.bed
wc -l dados/bedtools/AMOSTRA01_merged_sorted.bed
```
_bedtools coverage para encontrar cobertura em cada região_

```bash
bedtools coverage -a dados/bedtools/AMOSTRA01_merged_sorted.bed \
-b dados/bwa/AMOSTRA01_sorted.bam -mean \
>dados/bedtools/AMOSTRA01_coverageBed_mean.bed

head dados/bedtools/AMOSTRA01_coverageBed_mean.bed

tree dados
```
```jupyter
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

tabela = pd.read_csv("dados/bedtools/AMOSTRA01_coverageBed_mean.bed",header=None, names=np.array(["chr","start","end","cobertura"]), sep="\t")
tabela["cobertura"].describe()
```
##### CHAMADA DE VARIANTES COM O GATK

```bash
time gatk-4.2.2.0/gatk HaplotypeCaller \
-R referencia/grch38/grch38.chr22.fasta \
-I dados/bwa/AMOSTRA01_sorted.bam \
-O dados/gatk/AMOSTRA01_sorted.vcf \
-bamout dados/bwa/AMOSTRA01_sorted_bamout.bam

tree dados

cat dados/gatk/AMOSTRA01_sorted.vcf | grep -v "'#" | wc -l
```
##### ANOTAÇÃO DAS VARIANTES