# Genome Mapping Experiment of Snail/Twist KO to LV, trimmed with TrimGalore and mapped with STAR

1. Index Genome Reference

```bash
# Go to the Interactive Node:

srun -p interactive --pty bash
```
Indexing Lithechinus variegatus genome

```bash
srun -p interactive --pty /bin/bash
mkdir STAR_genome
nano indexing.sh
#!/usr/bin/env bash
#SBATCH --mem 50000
module load STAR
STAR --runMode genomeGenerate --runThreadN 16  --genomeDir STAR_Genome --genomeFastaFiles Lvar_scaffolds.fasta --sjdbGTFfile Lvar.braker.pasa.gff --genomeSAindexNbases 13


sbatch indexing.sh
```


2. Run fasta quality check with fastqc


```bash

module load fastqc

ls *.fastq.gz > reads.list
for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.fastqc.sh;
echo "fastqc $i" >> $root.fastqc.sh
done

for file in *fastqc.sh ; do sbatch $file ; done
```
3. By looking at the results of fastqc, I found illumina universal adapters were used. So, I will remove them with trimgalore.

```bash

module load cutadapt
module load TrimGalore/0.6.5-fasrc01

mkdir trim_galore

ls *R1_001.fastq.gz > readsR1.list
for i in `cat readsR1.list`; do
root=`basename $i _R1_001.fastq.gz`;
echo '#!/usr/bin/env bash' > $root.trimgalore.sh;
echo "trim_galore --illumina --paired --fastqc -o trim_galore/ ${root}_R1_001.fastq.gz ${root}_R2_001.fastq.gz " >> $root.trimgalore.sh
done


for file in *trimgalore.sh ; do sbatch $file ; done
```
