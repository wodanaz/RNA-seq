
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
