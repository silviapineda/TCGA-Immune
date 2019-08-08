import glob

###Read the file in gtf 
filenames = glob.glob("/Users/Pinedasans/TCGA-Immune/Data/Pancreas_Validation/report/*alignments_report.txt")
print("start")
fout = open("/Users/Pinedasans/TCGA-Immune/Data/Pancreas_Validation/total_reads.txt",'w')
for f in filenames:
    print(f)
    sample = f[59:95]
    print(sample)
    fh = open(f)
    for line in fh:
        line=line.rstrip()
        if line.startswith("Total sequencing reads"):
            words = line.split()
            print(words[3])
            fout.write(sample + ';' + words[3] + '\n')
            continue
print("Done")