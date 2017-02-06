#! /bin/sh

# gm=../../graphmap/bin/Linux-x64/graphmap
gm=../../graphmap/bin/graphmap-not_release
gmspliced=../../graphmap-spliced/bin/graphmap-not_release

### Coverage 20x
# ref=reference-genomes/scerevisiae_S288c_with_mito_uppercase.fna
# gtf=reference-genomes/scerevisiae.gtf
# reads=reads-simulated/PacBio-cov20/transcriptome_scerevisiae_graphmap_minlen100/reads.fq 
# folder=evaluation/reads-simulated/PacBio-cov20/transcriptome_scerevisiae_graphmap_minlen100

### 10k reads
ref=reference-genomes/scerevisiae_S288c_with_mito_uppercase.fna
gtf=reference-genomes/scerevisiae.gtf
reads=reads-simulated/PacBio-10k/transcriptome_scerevisiae_graphmap_minlen100/reads.fq
folder=evaluation/reads-simulated/PacBio-10k/transcriptome_scerevisiae_graphmap_minlen100
mkdir -p $folder

### Run GraphMap on erroneous simulated reads with transcriptome mapping.
out=$folder/graphmap-trans-rna-v1
# /usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${out}.memtime \
echo "	$gm align -r $ref --gtf $gtf -d $reads -o ${out}.sam --extcigar 2>&1 | tee ${out}.tee"

### Run GraphMap on erroneous simulated reads with spliced approach.
out=$folder/graphmap-spliced-rna-v1
# /usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${out}.memtime \
# 	$gmspliced align --approach spliced -r $ref -d $reads -o ${out}.sam --extcigar 2>&1 | tee ${out}.tee

### Evaluate the results.
# ./run-evaluation-trans.py rna-v1

### Run GraphMap on error-corrected simulated reads.
### Test on error corrected reads.
# paf=reads-simulated/PacBio-cov20/transcriptome_scerevisiae_graphmap_minlen100/erc-reads.paf
# erc_reads=reads-simulated/PacBio-cov20/transcriptome_scerevisiae_graphmap_minlen100/erc-reads.fasta
# tools/minimap/minimap -Sw5 -L100 -m0 $reads $reads > $paf
# /home/isovic/work/eclipse-workspace/git/racon/bin/racon --erc $reads $paf $reads ${erc_reads}
# out=/home/isovic/work/eclipse-workspace/git/aligneval/evaluation/reads-simulated/PacBio-cov20/transcriptome_scerevisiae_graphmap_minlen100/graphmap-erc-v1.sam
# $gm align -r $ref --gtf $gtf -d ${erc_reads} -o $out
