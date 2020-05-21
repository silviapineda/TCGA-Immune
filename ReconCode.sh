#!/bin/bash

FILES=/Users/Pinedasans/TCGA-Immune/Data/Validation_Normal_pancreas/Recon/*.txt
for f in $FILES
do  
  echo $f
  python2.7 ~/programs/Recon-master2/recon_v2.2.py -R -c -t 30 -o $f.fitfile.txt $f
done


FILES=/Users/Pinedasans/TCGA-Immune/Data/Validation_Normal_pancreas/Recon/*fitfile.txt

variantParameters=
for f in $FILES
do
    echo $f
    variantParameters="$variantParameters $f"
done

echo $variantParameters
python2.7 ~/programs/Recon-master2/recon_v2.2.py -D -Q 0 1 2 inf -b ~/programs/Recon-master2/error_bar_parameters.txt -o /Users/Pinedasans/TCGA-Immune/Data/Validation_Normal_pancreas/Recon/text_D_number.txt $variantParameters

