## 1. Align sequencing reads
FILES="list of SRR"
for i in $FILES]
do
mixcr align -p rna-seq -s hsa -f -r SRR1069352_aligments_report.txt -t 20 -OallowPartialAlignments=true SRR1069352_pass_1.fastq.gz SRR1069352_pass_2.fastq.gz SRR1069352_alignments.vdjca
done

#2. Partial alignments
for i in $(ls *alignments.vdjca | rev | cut -c 17- | rev | uniq)

do
echo $i 

#mixcr assemblePartial -f -r SRR1069352_aligments_report_rescued_1.txt SRR1069352_alignments.vdjca SRR1069352_alignments_rescued_1.vdjca
#mixcr assemblePartial -f -r SRR1069352_aligments_report_rescued_2.txt SRR1069352_alignments_rescued_1.vdjca SRR1069352_alignments_rescued_2.vdjca

mixcr assemblePartial -f -r ${i}aligments_report_rescued_1.txt ${i}alignments.vdjca ${i}alignments_rescued_1.vdjca
mixcr assemblePartial -f -r ${i}aligments_report_rescued_2.txt ${i}alignments_rescued_1.vdjca ${i}alignments_rescued_2.vdjca


done

##3. Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr extend -f -r ${i}aligments_report_extended.txt  ${i}alignments_rescued_2.vdjca ${i}alignments_extended.vdjca
#mixcr extend -f -r SRR1069352_aligments_report_extended.txt  SRR1069352_alignments_rescued_2.vdjca SRR1069352_alignments_extended.vdjca

done

##4. Assemble the clonotypes
#for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

#do
#echo $i
#mixcr assemble -r ${i}clonotypes_report.txt -i ${i}index_file ${i}alignments_extended.vdjca ${i}output.clns

#done

##5. Extract Alignments
for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)

do
echo $i
mixcr exportAlignments -readId -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene ${i}alignments_extended.vdjca ${i}alignments.txt
#mixcr exportAlignments -f -readId -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene SRR1069352_alignments_extended.vdjca SRR1069352_alignments.txt

done

##6. Extract Clones
#for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)
#do
#echo $i
#mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${i}output.clns ${i}clones.txt

#mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3  {SRR3478950}_output.clns {SRR3478950}_clones.txt
#done