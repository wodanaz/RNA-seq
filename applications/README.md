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
STAR --runMode genomeGenerate --runThreadN 16  --genomeDir STAR_genome --genomeFastaFiles L_var_genome.fa --sjdbGTFfile L_var_noseq.gff --genomeSAindexNbases 13


sbatch indexing.sh
```

```


2. Run fasta quality check with fastqc

```
module load fastqc



3. Remove Nextera Adapters


```bash
ls *.fastq.gz > reads.list
for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.fastqc.sh;
echo "fastqc $i" >> $root.fastqc.sh
done

for file in *fastqc.sh ; do sbatch $file ; done
```

```bash

module load cutadapt
module load TrimGalore/0.6.5-fasrc01



ls *.fastq.gz > reads.list
for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.trimgalore.sh;
echo "trim_galore --fastqc --nextera $i  " >> $root.trimgalore.sh
done


trimmomatic PE LW-S01_1C1_S74_L005_R1_001.fastq.gz LW-S01_1C1_S74_L005_R2_001.fastq.gz trimmed_1.fq unpaired_1.fq trimmed_2.fq unpaired_2.fq SLIDINGWINDOW:4:30 TRAILING:30 ILLUMINACLIP:adapter.fa:2:30:5 


for file in *trimgalore.sh ; do sbatch $file ; done
```
