# Genome Mapping Experiment of Snail/Twist KO to LV, trimmed with TrimGalore and mapped with STAR

1. Index Genome Reference

```bash
# Go to the Interactive Node:

srun -p interactive --pty bash
```
Indexing Lithechinus variegatus genome

```bash


module load Anaconda3
conda activate /data/wraycompute/alejo/conda/aleconda


mkdir STAR_Genome


nano indexing.sh
#!/usr/bin/env bash
#SBATCH --mem 50000
module load STAR
STAR --runMode genomeGenerate --runThreadN 16  --genomeDir STAR_Genome --genomeFastaFiles STAR_Genome/Lvar_scaffolds.fasta --sjdbGTFfile STAR_Genome/Lvar.braker.pasa.gff --genomeSAindexNbases 13 --sjdbOverhang=150


sbatch indexing.sh
```


2. Run fasta quality check with fastqc


```bash

module load Anaconda3
conda activate /data/wraycompute/alejo/conda/aleconda



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

module load Anaconda3
conda activate /data/wraycompute/alejo/conda/aleconda



mkdir trim_galore

ls *R1_001.fastq.gz > readsR1.list
for i in `cat readsR1.list`; do
root=`basename $i _R1_001.fastq.gz`;
echo '#!/usr/bin/env bash' > $root.trimgalore.sh;
echo "trim_galore --illumina --paired --fastqc -o trim_galore/ ${root}_R1_001.fastq.gz ${root}_R2_001.fastq.gz " >> $root.trimgalore.sh
done


for file in *trimgalore.sh ; do sbatch $file ; done
```




4. Map trimmed reads using STAR


```bash

module load Anaconda3
conda activate /data/wraycompute/alejo/conda/aleconda



ls trim_galore/*_val_1.fq.gz > reads.list

for i in `cat reads.list`; do
root=`basename $i _R1_001_val_1.fq.gz`;
echo '#!/usr/bin/env bash' > $root.trim2bam.sh;
echo "#SBATCH --mail-type=END" >> $root.trim2bam.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> $root.trim2bam.sh;
echo "#SBATCH -c 16" >> $root.trim2bam.sh;
echo "#SBATCH --mem 20000" >> $root.trim2bam.sh;
echo "module load STAR" >> $root.trim2bam.sh;
echo "STAR --runThreadN 16 --outFilterMismatchNoverLmax 0.05 --genomeDir /data/wraycompute/alejo/bulk_RNAseq/genome/STAR_Genome --readFilesIn $i trim_galore/${root}_R2_001_val_2.fq.gz --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --twopassMode Basic --outReadsUnmapped Fastx --outFileNamePrefix ${root}." >> $root.trim2bam.sh
done

for file in LW*trim2bam.sh ; do sbatch $file ; done



```




4. Obtain counts using HTSeq and samtools

```bash
module load HTSeq
module load samtools

ls *.bam > bams.list
for i in `cat bams.list`; do
root=`basename $i .Aligned.sortedByCoord.out.bam`;
echo '#!/usr/bin/env bash' > $root.bam2count.sh;
echo "#SBATCH --mail-type=END" >> $root.bam2count.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> $root.bam2count.sh;
echo "#SBATCH --mem 15000 " >> $root.bam2count.sh;
echo "module load HTSeq" >> $root.bam2count.sh;
echo "htseq-count --format=bam --stranded=no --type=gene --order=pos --idattr=ID $i STAR_Genome/Lvar.braker.pasa.gff > $root.counts.txt" >> $root.bam2count.sh
done

for file in *bam2count.sh ; do sbatch $file ; done


```

5. Compile data into a single table


```bash

perl compiler.pl *.counts.txt > allcounts.txt

cp allcounts.txt allcounts.backup.txt

grep "__" allcounts.txt -v > snail_twist_counts.txt



sed -r 's/ /\t/g' snail_twist_counts.txt | sed -r 's/LW-//g' | sed -r 's/.counts.txt//g' > snail.twist_counts_final.txt



```


