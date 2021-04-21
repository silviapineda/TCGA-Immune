## 1. Align sequencing reads

#!/bin/bash/
echo "start"
for i in $(ls FASTQ/*fastq.gz | rev | cut -c 12- | rev | uniq)
do
f=$(basename "$i")
echo $f
/home/spineda/miniconda3/bin/mixcr align -p rna-seq -s hsa -f -r ${f}_alignments_report.txt -t 20 -OallowPartialAlignments=true ${i}_1.
fastq.gz ${i}_2.fastq.gz ${f}_alignments.vdjca
done
echo "end"

##2. Partial alignments

#!/bin/bash/

echo "start"
for i in $(ls FASTQ/*alignments.vdjca | rev | cut -c 18- | rev | uniq)
do
f=$(basename "$i")
echo $f
/home/spineda/miniconda3/bin/mixcr assemblePartial -f -r ${f}_alignments_report_rescued_1.txt ${f}_alignments.vdjca ${f}_alignments_res
cued_1.vdjca
/home/spineda/miniconda3/bin/mixcr assemblePartial -f -r ${f}_alignments_report_rescued_2.txt ${f}_alignments_rescued_1.vdjca ${f}_alig
nments_rescued_2.vdjca
done
echo "end"


##3. Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences

#!/bin/bash/

echo "start"
for i in $(ls FASTQ/*alignments_rescued_2.vdjca | rev | cut -c 28- | rev | uniq)
do
f=$(basename "$i")
echo $f
/home/spineda/miniconda3/bin/mixcr extend -f -r ${f}_alignments_report_extended.txt  ${f}_alignments_rescued_2.vdjca ${f}_alignments_ex
tended.vdjca
done
echo "end"


##4. Extract Alignments

#!/bin/bash/

echo "start"
for i in $(ls FASTQ/*alignments_extended.vdjca | rev | cut -c 27- | rev | uniq)
do
f=$(basename "$i")
echo $f
/home/spineda/miniconda3/bin/mixcr exportAlignments -f -targetSequences -readIds -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dA
lignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene ${f}_alignments_extended.vdjca ${f}_alignments.txt
done
echo "end"