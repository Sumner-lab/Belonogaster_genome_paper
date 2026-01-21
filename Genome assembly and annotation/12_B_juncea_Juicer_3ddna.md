# Scaffolding B. juncea genome with juicer and Hi-C

### extract reads
```
/work/gif/archiveNova/2022_TothBeeAssemblies/IOW2052_B_juncea
#use huge node
module load picard/2.17.0-ft5qztz;  java  -Xmx900G -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar SamToFastq I=DTG-OmniC-219_R1_001.fastq.gz.bam  FASTQ=DTG-OmniC-219_R1_001_R1.fq F2=DTG-OmniC-219_R1_001_R2.fq FU=DTG-OmniC-219_R1_001__unpaired.fq
```

### Setup juicer
```
/work/gif/remkv6/Toth/07_B_juncea/01_Juicer
cp -rf /work/gif/00_JuicerSkeleton .
mv 00_JuicerSkeleton 01_Juicer

/work/gif/remkv6/Toth/07_B_juncea/01_Juicer/references
#rename scaffolds so they are more managable
awk '/^>/{print ">Scaffold_" ++i; next}{print}'  /work/gif/archiveNova/2022_TothBeeAssemblies/IOW2052_B_juncea/jasmine-iow2052-mb-hirise-ju840__01-07-2022__hic_output.fasta >B_junceaGenome.fasta
cd ../
bioawk -c fastx '{print $name"\t"length($seq)}' references/B_junceaGenome.fasta >chrom.sizes
cd fastq/
ln -s /work/gif/archiveNova/2022_TothBeeAssemblies/IOW2052_B_juncea/DTG-OmniC-219_R1_001_R1.fq Bjuncea_R1.fastq
ln -s /work/gif/archiveNova/2022_TothBeeAssemblies/IOW2052_B_juncea/DTG-OmniC-219_R1_001_R2.fq Bjuncea_R2.fastq
cd ../splits
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/Bjuncea_R1.fastq &
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/Bjuncea_R2.fastq &

/work/gif/remkv6/Toth/06_L_flavolineata/02_Juicer

echo "ml jdk; ml bwa; ml juicer/1.5.7; juicer.sh -d /work/gif/remkv6/Toth/07_B_juncea/01_Juicer -p /work/gif/remkv6/Toth/07_B_juncea/01_Juicer/chrom.sizes -s none -z /work/gif/remkv6/Toth/07_B_juncea/01_Juicer/references/B_junceaGenome.fasta  -q short -Q 2:00:00 -l medium -L 12:00:00 -t 36" >juicer.sh
```

### Run 3ddna
```
/work/gif/remkv6/Toth/07_B_juncea/01_Juicer/3Ddna/01_3d-dnaDefault
ln -s ../../aligned/merged_nodups.txt
ln -s ../../references/B_junceaGenome.fasta
module load miniconda3/4.3.30-qdauveb;source activate 3d-dna;module load jdk;module load parallel;cd /work/gif/remkv6/Toth/07_B_juncea/01_Juicer/3Ddna/01_3d-dnaDefault ;bash run-asm-pipeline.sh /work/gif/remkv6/Toth/07_B_juncea/01_Juicer/3Ddna/01_3d-dnaDefault/B_junceaGenome.fasta /work/gif/remkv6/Toth/07_B_juncea/01_Juicer/3Ddna/01_3d-dnaDefault/merged_nodups.txt

```

### Assembly statistics
Received Assembly
```
---------------- Information for assembly 'B_junceaGenome.fasta' ----------------


                                         Number of scaffolds        206
                                     Total size of scaffolds  315721248
                                            Longest scaffold   11540610
                                           Shortest scaffold       4162
                                 Number of scaffolds > 1K nt        206 100.0%
                                Number of scaffolds > 10K nt        204  99.0%
                               Number of scaffolds > 100K nt        157  76.2%
                                 Number of scaffolds > 1M nt         65  31.6%
                                Number of scaffolds > 10M nt          4   1.9%
                                          Mean scaffold size    1532627
                                        Median scaffold size     382434
                                         N50 scaffold length    4927384
                                          L50 scaffold count         21
                                         n90 scaffold length     835833
                                          L90 scaffold count         72
                                                 scaffold %A      33.93
                                                 scaffold %C      16.06
                                                 scaffold %G      15.99
                                                 scaffold %T      34.02
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      51.1%
              Percentage of assembly in unscaffolded contigs      48.9%
                      Average number of contigs per scaffold        1.2
Average length of break (>25 Ns) between contigs in scaffold        100

                                           Number of contigs        237
                              Number of contigs in scaffolds         54
                          Number of contigs not in scaffolds        183
                                       Total size of contigs  315718148
                                              Longest contig    7306445
                                             Shortest contig       4162
                                   Number of contigs > 1K nt        237 100.0%
                                  Number of contigs > 10K nt        235  99.2%
                                 Number of contigs > 100K nt        183  77.2%
                                   Number of contigs > 1M nt         86  36.3%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size    1332144
                                          Median contig size     448708
                                           N50 contig length    3793124
                                            L50 contig count         32
                                           n90 contig length     790695
                                            L90 contig count         98
                                                   contig %A      33.93
                                                   contig %C      16.06
                                                   contig %G      15.99
                                                   contig %T      34.02
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```
scaffolded assembly
```
---------------- Information for assembly 'B_junceaGenome.FINAL.fasta' ----------------


                                         Number of scaffolds         99
                                     Total size of scaffolds  315790748
                                            Longest scaffold   25458634
                                           Shortest scaffold       1000
                                 Number of scaffolds > 1K nt         97  98.0%
                                Number of scaffolds > 10K nt         83  83.8%
                               Number of scaffolds > 100K nt         53  53.5%
                                 Number of scaffolds > 1M nt         39  39.4%
                                Number of scaffolds > 10M nt          8   8.1%
                                          Mean scaffold size    3189806
                                        Median scaffold size     133699
                                         N50 scaffold length    7911001
                                          L50 scaffold count         14
                                         n90 scaffold length    5357696
                                          L90 scaffold count         33
                                                 scaffold %A      34.26
                                                 scaffold %C      16.07
                                                 scaffold %G      15.97
                                                 scaffold %T      33.67
                                                 scaffold %N       0.02
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs      91.7%
              Percentage of assembly in unscaffolded contigs       8.3%
                      Average number of contigs per scaffold        2.7
Average length of break (>25 Ns) between contigs in scaffold        427

                                           Number of contigs        269
                              Number of contigs in scaffolds        206
                          Number of contigs not in scaffolds         63
                                       Total size of contigs  315718148
                                              Longest contig    7306445
                                             Shortest contig        616
                                   Number of contigs > 1K nt        265  98.5%
                                  Number of contigs > 10K nt        246  91.4%
                                 Number of contigs > 100K nt        188  69.9%
                                   Number of contigs > 1M nt         86  32.0%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size    1173673
                                          Median contig size     356350
                                           N50 contig length    3790535
                                            L50 contig count         32
                                           n90 contig length     703441
                                            L90 contig count        101
                                                   contig %A      34.27
                                                   contig %C      16.08
                                                   contig %G      15.97
                                                   contig %T      33.68
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

```


### Busco scores
Received assembly
```
/work/gif/remkv6/Toth/07_B_juncea/02_busco/01_receivedAssembly
ln -s ../../01_Juicer/references/B_junceaGenome.fasta
echo  "source activate busco5_env ; busco -i B_junceaGenome.fasta -o B_junceaGenome.Busco -m geno --auto-lineage-euk --long --augustus -c 35 -f"  >busco.sh

--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:100.0%[S:98.4%,D:1.6%],F:0.0%,M:0.0%,n:255     |
|255    Complete BUSCOs (C)                       |
|251    Complete and single-copy BUSCOs (S)       |
|4      Complete and duplicated BUSCOs (D)        |
|0      Fragmented BUSCOs (F)                     |
|0      Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset hymenoptera_odb10           |
--------------------------------------------------
|C:97.5%[S:97.1%,D:0.4%],F:0.6%,M:1.9%,n:5991     |
|5844   Complete BUSCOs (C)                       |
|5820   Complete and single-copy BUSCOs (S)       |
|24     Complete and duplicated BUSCOs (D)        |
|34     Fragmented BUSCOs (F)                     |
|113    Missing BUSCOs (M)                        |
|5991   Total BUSCO groups searched               |
--------------------------------------------------

```

Scaffolded Assembly
```
--------------------------------------------------
|Results from generic domain eukaryota_odb10      |
--------------------------------------------------
|C:100.0%[S:98.4%,D:1.6%],F:0.0%,M:0.0%,n:255     |
|255    Complete BUSCOs (C)                       |
|251    Complete and single-copy BUSCOs (S)       |
|4      Complete and duplicated BUSCOs (D)        |
|0      Fragmented BUSCOs (F)                     |
|0      Missing BUSCOs (M)                        |
|255    Total BUSCO groups searched               |
--------------------------------------------------

--------------------------------------------------
|Results from dataset hymenoptera_odb10           |
--------------------------------------------------
|C:97.4%[S:97.0%,D:0.4%],F:0.6%,M:2.0%,n:5991     |
|5838   Complete BUSCOs (C)                       |
|5814   Complete and single-copy BUSCOs (S)       |
|24     Complete and duplicated BUSCOs (D)        |
|34     Fragmented BUSCOs (F)                     |
|119    Missing BUSCOs (M)                        |
|5991   Total BUSCO groups searched               |
--------------------------------------------------

```