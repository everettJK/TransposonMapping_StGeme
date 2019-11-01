MariaDB [specimen_management]> select SpecimenAccNum, Patient, CellType, Timepoint, SpecimenInfo from gtsp where Trial = 'StGeme Tn';
+----------------+---------------+--------------------+-----------+----------------------------------------------------------------------+
| SpecimenAccNum | Patient       | CellType           | Timepoint | SpecimenInfo                                                         |
+----------------+---------------+--------------------+-----------+----------------------------------------------------------------------+
| GTSP3171       | pChallenge1.1 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 1, technical replicate 1" |
| GTSP3172       | pChallenge1.2 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 1, technical replicate 2" |
| GTSP3173       | pChallenge2.1 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 2, technical replicate 1" |
| GTSP3174       | pChallenge2.2 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 2, technical replicate 2" |
| GTSP3175       | pChallenge3.1 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 3, technical replicate 1" |
| GTSP3176       | pChallenge3.2 | Bacterial Cell Iso | min60     | "Challenge condition: biological replicate 3, technical replicate 2" |
| GTSP3177       | pControl1.1   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 1, technical replicate 1"   |
| GTSP3178       | pControl1.2   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 1, technical replicate 2"   |
| GTSP3179       | pControl2.1   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 2, technical replicate 1"   |
| GTSP3180       | pControl2.2   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 2, technical replicate 2"   |
| GTSP3181       | pControl3.1   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 3, technical replicate 1"   |
| GTSP3182       | pControl3.2   | Bacterial Cell Iso | min60     | "Control condition: biological replicate 3, technical replicate 2"   |
+----------------+---------------+--------------------+-----------+----------------------------------------------------------------------+


ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/470/375/GCF_000470375.1_Vir5453/GCF_000470375.1_Vir5453_genomic.fna.gz
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/470/375/GCF_000470375.1_Vir5453/GCF_000470375.1_Vir5453_genomic.gbff.gz

~/software/bwa/bwa index -p GCF_000470375.1_Vir5453 -a bwtsw ../GCF_000470375.1_Vir5453_genomic.fna
~/software/bwa/bwa mem -M -t 30 ../data/genome/BWA/GCF_000470375.1_Vir5453 ../data/Undetermined_S0_L001_R1_001.fastq.gz ../data/Undetermined_S0_L001_R2_001.fastq.gz > GCF_000470375.1_Vir5453.sam


https://www.ncbi.nlm.nih.gov/genome/3165?genome_assembly_id=259208



