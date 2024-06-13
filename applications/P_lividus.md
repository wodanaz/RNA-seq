
#To analyse new data from Paracentrotus, France

```bash
conda activate /data/wraycompute/alejo/conda/aleconda


ls *.fastq.gz > reads.list



for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.fastqc.sh;
echo "fastqc $i" >> $root.fastqc.sh
done



for file in *fastqc.sh ; do sbatch $file ; done


```



# For P. lividus
```bash
conda activate /data/wraycompute/alejo/conda/aleconda


for i in `cat reads.list`; do
root=`basename $i .fq.gz`;
echo '#!/usr/bin/env bash' > $root.trim2bam.sh;
echo "#SBATCH --mail-type=END" >> $root.trim2bam.sh;
echo "#SBATCH --mail-user=ab620@duke.edu" >> $root.trim2bam.sh; 
echo "#SBATCH -c 16" >> $root.trim2bam.sh;
echo "#SBATCH --mem 20000" >> $root.trim2bam.sh;
echo "STAR --runThreadN 16 --outFilterMismatchNoverLmax 0.1 --genomeDir /data/wraycompute/alejo/av-axis-Pl-v2/P_lividus/STAR_Genome  --readFilesIn $i --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --twopassMode Basic --outReadsUnmapped Fastx --outFileNamePrefix ${root}." >> $root.trim2bam.sh 
done

for file in *trim2bam.sh ; do sbatch $file ; done


```

# For H. erythrogramma


```bash
conda activate /data/wraycompute/alejo/conda/aleconda


for i in `cat reads.list`; do
root=`basename $i .fq.gz`;
echo '#!/usr/bin/env bash' > $root.trim2bam.sh;
echo "#SBATCH --mail-type=END" >> $root.trim2bam.sh;
echo "#SBATCH --mail-user=ab620@duke.edu" >> $root.trim2bam.sh; 
echo "#SBATCH -c 16" >> $root.trim2bam.sh;
echo "#SBATCH --mem 20000" >> $root.trim2bam.sh;
echo "STAR --runThreadN 16 --outFilterMismatchNoverLmax 0.1 --genomeDir /data/wraycompute/alejo/newgenomes/P_lividus/STAR_Genome  --readFilesIn $i --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --twopassMode Basic --outReadsUnmapped Fastx --outFileNamePrefix ${root}." >> $root.trim2bam.sh 
done

for file in *trim2bam.sh ; do sbatch $file ; done


```




