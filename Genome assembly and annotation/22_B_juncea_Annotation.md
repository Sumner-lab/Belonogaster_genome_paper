# Annotate the B. juncea genome


### RNA-seq alignment
```
/work/gif/remkv6/Toth/07_B_juncea/04_RNAseqAlignment

for f in  /work/gif/archiveNova/2022_TothBeeAssemblies/IOW2052_B_juncea/RNA-seq/raw_data/*/*gz ;do ln -s $f;done
ln -s ../01_Juicer/3Ddna/04_3d-dnaDiploidRepCovFinalize/B_junceaGenome.FINAL.fasta

ml star;STAR  --runMode genomeGenerate --genomeDir /work/gif/remkv6/Toth/07_B_juncea/04_RNAseqAlignment/ --genomeFastaFiles B_junceaGenome.FINAL.fasta

for f in *_1.fq.gz; do echo "ml star;STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --genomeDir /work/gif/remkv6/Toth/07_B_juncea/04_RNAseqAlignment/ --outFileNamePrefix  "${f%*_1.fq.gz}" --readFilesIn  "$f" "${f%*_1.fq.gz}"_2.fq.gz";done >AlignStar.sh

# merge the bams
ls -lrth *out.bam |awk '$5!="0"{print $9}' >bam.list
ml samtools; samtools merge -@ 24  -b bam.list MergedRNA.bam ;samtools sort -o MergedRNA_sorted.bam -T TEMP  --threads 24   MergedRNA.bam
```

### Call repeats in the genome
```
/work/gif/remkv6/Toth/07_B_juncea/05_Repeats
 ln -s ../02_busco/02_ScaffoldedAssembly/B_junceaGenome.FINAL.fasta


ml miniconda3; source activate repeatmodeler2;
BuildDatabase -name coarctatus B_junceaGenome.FINAL.fasta
ml miniconda3; source activate repeatmodeler2; RepeatModeler -database coarctatus -pa 36
ln -s */consensi.fa.classified
ml miniconda3; source activate repeatmasker; RepeatMasker -pa 36 -norna  -dir RepeatmaskerOut -small -gff -lib consensi.fa.classified B_junceaGenome.FINAL.fasta
#softmask
ml bedtools2;bedtools maskfasta -fi ../02_busco/02_ScaffoldedAssembly/B_junceaGenome.FINAL.fasta -fo SoftmaskedB_junceaGenome.FINAL.fasta -bed RepeatmaskerOut/B_junceaGenome.FINAL.fasta.out.gff  -soft
```


### Run Braker
```
/work/gif/remkv6/Toth/07_B_juncea/06_Braker
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus;
cp -rf /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/augustus-3.3.2-rtcnsefyulxnfscwgpwi5tc7civqwvsq/bin/ .
cd scripts
cp -rf  /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/braker-2.1.2-75wblifp2zieps5rf7tzp7ajcwvzo2oz/cfg/ .
cd ../config/species/
cp -rf ../../../../02_busco/02_ScaffoldedAssembly/B_junceaGenome.FINAL.Busco/run_hymenoptera_odb10/augustus_output/retraining_parameters/BUSCO_B_junceaGenome.FINAL.Busco .

ln -s ../05_Repeats/SoftmaskedB_junceaGenome.FINAL.fasta


echo "ml samtools; samtools index MergedRNA_sorted.bam;ml bamtools;ml genemark-et/4.38-63ipkx4;ml augustus/3.3.2-py3-openmpi3-gee2bjt;ml genemark-et/4.38-63ipkx4; ml braker/2.1.2-py3-openmpi3-75wblif;braker.pl --genome=SoftmaskedB_junceaGenome.FINAL.fasta --softmasking --species=BUSCO_B_junceaGenome.FINAL.Busco --bam=MergedRNA_sorted.bam  --AUGUSTUS_CONFIG_PATH=/work/gif/remkv6/Toth/07_B_juncea/06_Braker/Augustus/config/ --overwrite --useexisting" >juncBraker.sh

ml braker;ml augustus;cat augustus.hints.gtf |gtf2gff.pl --gff3 --out=B_junceaBraker.gff3

```

Annotation Stats
```
awk '$3=="gene"{print $5-$4}' B_junceaBraker.gff3|summary.sh
Total:  100,282,108
Count:  18,512
Mean:   5,417
Median: 1,683
Min:    123
Max:    338,422
awk '$3=="mRNA"{print $5-$4}' B_junceaBraker.gff3|summary.sh
Total:  140,866,210
Count:  20,723
Mean:   6,797
Median: 1,940
Min:    123
Max:    338,422
awk '$3=="CDS"{print $5-$4}' B_junceaBraker.gff3|summary.sh
Total:  31,939,949
Count:  124,821
Mean:   255
Median: 179
Min:    2
Max:    13,756
```


### Busco on Annotation
```
/work/gif/remkv6/Toth/07_B_juncea/02_busco/03_BrakerAnnotation
ln -s ../../06_Braker/braker/augustus.hints.aa
sed 's/\*//g' augustus.hints.aa   >juncBraker_proteins.fasta

sh runBusco.sh juncBraker_proteins.fasta

#runBusco
#############################################################
#!/bin/bash
#runBusco.sh
#Here is how to run this script: sh runBusco.sh Genome.fasta
Genome="$1"

ml miniconda3;source activate busco5_env ;
busco -i ${Genome} \
-o ${Genome%.*}_Busco \
-m prot \
--auto-lineage-euk \
-c 35 \
-f
##############################################################
```

Results
```
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:98.4%[S:90.2%,D:8.2%],F:0.8%,M:0.8%,n:255      |
|251    Complete BUSCOs (C)                       |
|230    Complete and single-copy BUSCOs (S)       |
|21     Complete and duplicated BUSCOs (D)        |
|2      Fragmented BUSCOs (F)                     |
|2      Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset hymenoptera_odb10           |
--------------------------------------------------
|C:95.3%[S:79.0%,D:16.3%],F:0.9%,M:3.8%,n:5991    |
|5708   Complete BUSCOs (C)                       |
|4734   Complete and single-copy BUSCOs (S)       |
|974    Complete and duplicated BUSCOs (D)        |
|52     Fragmented BUSCOs (F)                     |
|231    Missing BUSCOs (M)                        |
|5991   Total BUSCO groups searched               |
--------------------------------------------------

```