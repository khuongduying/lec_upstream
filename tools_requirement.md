# INSTALL TOOLS

## List of tools 
| Tools | Description | file format | Preferences |
|--- | ------------- | ---------------- | -------------- |
| FastQC | Sequencing quality control | fastq | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Trimmomatic | Useful trimming tasks for illumina paired-end and single ended data. | fastq | https://github.com/usadellab/Trimmomatic |
| BWA mem | Mapped reads to reference genome | sam/bam | https://bio-bwa.sourceforge.net/bwa.shtml |
| picard | | |
| IGV |  Vsual exploration of genomic data | cram | https://igv.org/
| GATK | A genomic analysis toolkit focused on variant discovery. | | https://gatk.broadinstitute.org/hc/en-us |
| samtools | Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format | sam/bam/cram | http://www.htslib.org/ |

## FastQC
Download: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
```bash
#Install & unzip
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip

unzip fastqc_v0.12.1.zip

# Write at the end of the .bashrc file. This command will let you export the path to FastQC, and execute it everywhere
nano ~/.bashrc

export PATH='path/to/FastQC/':$PATH
# For example: 
export PATH='/home/duydao/dnaseq_work/tools/FastQC':$PATH

source ~/.bashrc

#Try it by running
fastqc
```

## Trimmomatic
Download: https://anaconda.org/bioconda/trimmomatic/0.39/download/noarch/trimmomatic-0.39-hdfd78af_2.tar.bz2

```bash
# Move the file to tools/ and then unzip
mkdir trimmomatic
tar -xvf trimmomatic-0.39-hdfd78af_2.tar.bz2 -C trimmomatic/
#
nano ~/.bashrc
#
export PATH='/home/duydao/dnaseq_work/tools/trimmomatic/bin/':$PATH
#
source ~/.bashrc
```
## BWA
```bash
git clone https://github.com/lh3/bwa.git

cd bwa/

make

# export to PATH
nano ~/.bashrc

export PATH='/home/duydao/dnaseq_work/tools/bwa/':$PATH #Write at the end of the file.

source ~/.bashrc

# done.
```

## Picard
Download: https://github.com/broadinstitute/picard/zipball/master
```bash
# Move picard to tools/ and then unzip
unzip
cd 
./gradlew bundle

```
https://github.com/broadinstitute/picard/blob/master/README.md#building-picard

## GATK
https://github.com/broadinstitute/gatk/releases

```bash
git clone https://github.com/broadinstitute/gatk.git

# Install git-lfs
wget https://github.com/git-lfs/git-lfs/releases/download/v3.3.0/git-lfs-linux-amd64-v3.3.0.tar.gz

tar -xvf git-lfs-linux-amd64-v3.3.0.tar.gz

# Set up java version 17
## Download java 17
sudo apt install openjdk-17-jdk

## Switch to java 17
sudo update-alternatives --config java
## Type selection number 1,2,3,... that correspond to Java 17 and press Enter.
  Selection    Path                                         Priority   Status
------------------------------------------------------------
  0            /usr/lib/jvm/java-19-openjdk-amd64/bin/java   1911      auto mode
* 1            /usr/lib/jvm/java-17-openjdk-amd64/bin/java   1711      manual mode
  2            /usr/lib/jvm/java-18-openjdk-amd64/bin/java   1811      manual mode
  3            /usr/lib/jvm/java-19-openjdk-amd64/bin/java   1911      manual mode

--> 1


# Build gatk
cd gatk/
sudo ./gradlew bundle # This may take a while.

# export to PATH
nano ~/.bashrc
#
export PATH='/home/duydao/dnaseq_work/tools/gatk/':$PATH #Write at the end of the file.
#
source ~/.bashrc
# done.
```

## Samtools
```bash
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2

tar -xvf samtools-1.17.tar.bz2 

cd samtools-1.17/

./configure

make

make install

# export to PATH
nano ~/.bashrc
#
export PATH='/home/duydao/dnaseq_work/tools/samtools-1.17/':$PATH #Write at the end of the file.
#
source ~/.bashrc
# done.
```



## Tabix
```bash
wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2

tar -xvf htslib-1.17.tar.bz2

cd htslib-1.17/
make

# export to PATH
nano ~/.bashrc
#
export PATH='/home/duydao/dnaseq_work/tools/htslib-1.17/':$PATH #Write at the end of the file.
#
source ~/.bashrc
# done.
```

## When complete installing all the tools, our tools dir will look like this
```bash
tree -L 1 tools

#
tools
├── bwa
├── FastQC
├── gatk
├── git-lfs-3.3.0
├── htslib-1.17
├── samtools-1.17
└── trimmomatic
```