
##Download the data
##First, we need to configure the sratoolkit 
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=dbgap_use


##prefetch
/opt/bioinfo/sratoolkit/bin/prefetch --list cart_DAR81372_201905280938.krt #This is to show what is in the cart file with all the ID's

FILES=$( cat srr_list.txt)

for i in $FILES
do
##To download without a cart you need to be in the /dbGap-21956/sra/ folder
/opt/bioinfo/sratoolkit/bin/prefetch --max-size 30G ${i} ## --max size is to increase the limit that can be downloaded at once, otherwise the files cannot be download

##fastq-dump
#https://edwards.sdsu.edu/research/fastq-dump/
/opt/bioinfo/sratoolkit/bin/fastq-dump --gzip --split-files --readids -F --skip-technical -W --read-filter pass ${i}.sra
done

#### Material and Methods from GTEX data
#https://science.sciencemag.org/content/sci/suppl/2015/05/06/348.6235.648.DC1/GTEx.SM.pdf

FILES="list of SRR"

## 1. Align sequencing reads
for i in $FILES
do
mixcr align -p rna-seq -s hsa -f -r ${i}_aligments_report.txt -t 20 -OallowPartialAlignments=true ${i}_pass_1.fastq.gz ${i}_pass_2.fastq.gz ${i}_alignments.vdjca
done

#2. Partial alignments
for i in $(ls *alignments.vdjca | rev | cut -c 17- | rev | uniq)

do
echo $i 

mixcr assemblePartial -f -r ${i}aligments_report_rescued_1.txt ${i}alignments.vdjca ${i}alignments_rescued_1.vdjca
mixcr assemblePartial -f -r ${i}aligments_report_rescued_2.txt ${i}alignments_rescued_1.vdjca ${i}alignments_rescued_2.vdjca

done

##3. Perform extension of incomplete TCR CDR3s with uniquely determined V and J genes using germline sequences
for i in $(ls *alignments_rescued_2.vdjca | rev | cut -c 27- | rev | uniq)

do
echo $i
mixcr extend -f -r ${i}aligments_report_extended.txt  ${i}alignments_rescued_2.vdjca ${i}alignments_extended.vdjca
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

done

##6. Extract Clones
#for i in $(ls *alignments_extended.vdjca | rev | cut -c 26- | rev | uniq)
#do
#echo $i
#mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${i}output.clns ${i}clones.txt

#mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3  {SRR3478950}_output.clns {SRR3478950}_clones.txt
#done