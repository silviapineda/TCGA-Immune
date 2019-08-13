#!/

bin/bash

FILES=/Users/Pinedasans/TCGA-Immune/Data/GTEx/Blood/RECON/*.txt
for f in $FILES
do  
   echo $f
   python2.7 ~/programs/Recon-master2/recon_v2.2.py -R -c -t 30 -o $f.fitfile.txt $f
   #python2.7 recon_v2.2.py -x --x_max 30 -o $f.plotfile.txt -b error_bar_parameters.txt $f.fitfile.txt
done


FILES=/Users/Pinedasans/TCGA-Immune/Data/GTEx/Blood/RECON/*fitfile.txt

variantParameters=
for f in $FILES
do
    echo $f
    variantParameters="$variantParameters $f"
done

echo $variantParameters
python2.7 ~/programs/Recon-master2/recon_v2.2.py -D -Q 0 1 2 inf -b error_bar_parameters.txt -o /:Users/Pinedasans/TCGA-Immune/Data/GTEx/Blood/RECON//test_D_number_table.txt $variantParameters

