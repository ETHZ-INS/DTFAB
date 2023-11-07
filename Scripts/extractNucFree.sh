#!/bin/bash
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
  echo "Usage: `basename $0` path/to/bam/folder path/to/output/folder [species]"
  echo "If species (mm/hs, or effective genome size) is given, peaks will be called on merged NF reads."
  exit
fi

mkdir -p $2

echo "Subsetting reads..."

for f in $1/*.bam; do
  echo $f
  f2=$2/`basename $f .bam`.NF.bam
  samtools view -h $f | \
  awk 'substr($0,1,1)=="@" || ($9>=30 && $9<=120) || ($9<=-30 && $9>=-120)' | \
  samtools view -b > $f2
  samtools index $f2
done

if [ -z "$3" ]; then
  echo "Done!"
  exit
else
  echo "Merging reads and calling peaks..."
  samtools merge -@ 4 $2/NFmerged.bam $2/*.NF.bam && \
    samtools index $2/NFmerged.bam
  macs2 callpeak -t $2/NFmerged.bam -n $2/merged.NF -f BAMPE -g $3 -q 0.01
fi

echo "Done!"
