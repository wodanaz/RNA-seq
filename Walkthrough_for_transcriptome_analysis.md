This Walkthrough generates count files using bowtie2 mapping tool from fastq files reconstructed from bam in order to run differential transcription analysis




```bash
ls -1 *bam > reads.list
for file in `cat reads.list`; do
root=`basename $file .sorted.bam`;
echo '#!/usr/bin/env bash' > bam2fq.$root.sh;
echo "#SBATCH --mail-type=END" >> bam2fq.$root.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> bam2fq.$root.sh;
echo "#SBATCH --mem 5000" >> bam2fq.$root.sh;
echo "module load samtools" >> bam2fq.$root.sh;
echo "samtools fastq $file | gzip > $root.fq.gz"  >> bam2fq.$root.sh;
done


for file in bam2fq*sh ; do sbatch $file ; done
```

################################
#------------------------------
# Quality check!!!! 

module load fastqc
mkdir fastqc_results
fastqc --outdir fastqc_results *fastq 



# error in 1625VP.sorted.bam and 5262HT.sorted.bam
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
[E::bgzf_read] Read block operation failed with error 4 after 0 of 4 bytes
[bam2fq_mainloop] Failed to read bam record.




# Download and unzip genome files
```bash
wget ftp://ftp.ensembl.org/pub/release-100/fasta/microtus_ochrogaster/dna/Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-100/gff3/microtus_ochrogaster/Microtus_ochrogaster.MicOch1.0.100.chr.gff3.gz
wget ftp://ftp.ensembl.org/pub/release-100/gtf/microtus_ochrogaster/Microtus_ochrogaster.MicOch1.0.100.chr.gtf.gz


gunzip Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa.gz
gunzip Microtus_ochrogaster.MicOch1.0.100.chr.gff3.gz
gunzip Microtus_ochrogaster.MicOch1.0.100.chr.gtf.gz


mkdir genome


cp Microtus_ochrogaster.MicOch1.0.100.chr.gtf genome
cp Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa genome
```


### INDEXING LV GENOME


```bash
srun -p interactive --pty /bin/bash
mkdir STAR_genome
nano indexing.sh
#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --mail-user=alebesc@gmail.com
#SBATCH -c 16
#SBATCH --mem 50000
module load STAR
STAR --runMode genomeGenerate --runThreadN 16  --genomeDir STAR_genome --genomeFastaFiles Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa --sjdbGTFfile Microtus_ochrogaster.MicOch1.0.100.chr.gtf


sbatch indexing.sh
```



```bash
module load STAR
ls -1 *fq.gz > fastq.list
for file in `cat fastq.list`; do
root=`basename $file .fq.gz`;
echo '#!/usr/bin/env bash' > Microtus.$root.sh;
echo "#SBATCH --mail-type=END" >> Microtus.$root.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> Microtus.$root.sh;
echo "#SBATCH -c 16" >> Microtus.$root.sh;
echo "#SBATCH --mem 30000" >> Microtus.$root.sh;
echo "module load STAR" >> Microtus.$root.sh;
echo "STAR --runThreadN 16 --outFilterMismatchNoverLmax 0.05 --genomeDir /data/wraycompute/alejo/phd2/genome/STAR_genome  --readFilesIn $file --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate  --readFilesCommand zcat  --twopassMode Basic --outReadsUnmapped Fastx --outFileNamePrefix ${root}." >> Microtus.$root.sh
done

for file in Microtus*sh ; do sbatch $file ; done
```


# 


```bash
module load HTSeq
module load samtools

ls *.bam > bams.list
for file in `cat bams.list`; do
root=`basename $file .Aligned.sortedByCoord.out.bam`;
echo '#!/usr/bin/env bash' > $root.bam2count.sh;
echo "#SBATCH --mail-type=END" >> $root.bam2count.sh;
echo "#SBATCH --mail-user=alebesc@gmail.com" >> $root.bam2count.sh;
echo "#SBATCH --mem 15000 " >> $root.bam2count.sh;
echo "module load HTSeq" >> $root.bam2count.sh;
echo "htseq-count --format=bam --stranded=no --type=gene --order=pos --idattr=gene_id $file Microtus_ochrogaster.MicOch1.0.100.chr.gtf > $root.counts.txt" >> $root.bam2count.sh
done



for file in *bam2count.sh ; do sbatch $file ; done
```



```bash
perl compiler.pl *.counts.txt > allcounts.txt

cp allcounts.txt allcounts.backup.txt




grep "__" allcounts.txt -v > vole_counts.txt



sed -r 's/ /\t/g' vole_counts.txt | sed -r 's/.counts.txt//g' > vole_counts_final.txt

```



# To Convert gtf to table to be used in R


```bash
 awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1] "\t"  a[3]  "\t" $1 ":" $4 "-" $5  }'  Microtus_ochrogaster.MicOch1.0.100.chr.gtf  |  sed 's/gene_id "//' | sed 's/"//g' | awk '{print $1 "\t" $3 "\t" $4 }' | sed -r 's/ensembl/none/g'  > moch_gene_annotations.tab

```
