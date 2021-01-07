# Genome Mapping Experiment of Snail/Twist KO to LV, trimmed TrimGalore and STAR

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
STAR --runMode genomeGenerate --runThreadN 16  --genomeDir STAR_genome --genomeFastaFiles Microtus_ochrogaster.MicOch1.0.dna.toplevel.fa --sjdbGTFfile Microtus_ochrogaster.MicOch1.0.100.chr.gtf


sbatch indexing.sh
```

```

2. Remove Nextera Adapters


```bash

module load cutadapt
module load TrimGalore/0.6.5-fasrc01



ls *.fastq.gz > reads.list
for i in `cat reads.list`; do
root=`basename $i .fastq.gz`;
echo '#!/usr/bin/env bash' > $root.trimgalore.sh;
echo "trim_galore --fastqc --nextera $i  " >> $root.trimgalore.sh
done


for file in *trimgalore.sh ; do sbatch $file ; done
```
