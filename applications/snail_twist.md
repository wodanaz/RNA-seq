# Assembling Viruses

### Genome Mapping Experiment of Snail/Twist KO to LV, trimmed TrimGalore and STAR

1. Index Genome Reference

```bash
# Go to the Interactive Node:

srun -p interactive --pty bash



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
