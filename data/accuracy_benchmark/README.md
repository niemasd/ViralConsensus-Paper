To benchmark the accuracy of ViralConsensus, simulated reads from known genomes and then compared the ViralConsensus sequence vs. the true genome sequence. We selected a high-confidence representative genome from each of the following SARS-CoV-2 lineages:
* **Lineage B:** [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)
* **Lineage B.1.1.7:** [LC650844.1](https://www.ncbi.nlm.nih.gov/nuccore/LC650844.1)
* **Lineage BA.1:** [OQ523614.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ523614.1)
* **Lineage BA.2:** [OQ194009.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ194009.1)
* **Lineage XBB.1:** [OQ346068.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ346068.1)

# Illumina
## Simulation
We used [ART version MountRainier-2016-06-05](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to simulate Illumina reads:

```bash
for f in lineage_* ; do for c in 30 40 50 ; do mkdir -p $f/c$c/illumina/fastq && for r in $(seq -w 1 10) ; do art_illumina -rs $RANDOM -q -na -ss HS20 -l 100 -f $c -i $f/*.fas -o $f/c$c/illumina/fastq/$f.c$c.illumina.r$r ; done ; done ; done
```

The individual ART command is as follows:

```bash
art_illumina -rs RNG_SEED -q -na -ss HS20 -l 100 -f COVERAGE -i REF_GENOME -o OUTPUT
```

* `-rs RNG_SEED` = Use `RNG_SEED` as the random number generation seed (`-rs $RANDOM` picks a random value from the OS)
* `-q` = Quiet mode
* `-na` = Don't output the alignment file (just the FASTQ)
* `-ss HS20` = Use the HiSeq 2000 error profile
* `-l 100` = Use a length of 100 bases for each read
* `-f COVERAGE` = Simulate at `COVERAGE` fold coverage (e.g. `-f 10` = 10X coverage)
* `-i REF_GENOME` = Use the reference genome in the FASTA file called `REF_GENOME`
* `-o OUTPUT` = Output file prefix

## Mapping
Illumina reads were then mapped to the reference genome using Minimap2's short-read preset and piped to Samtools to convert to BAM:

```bash
for f in lineage_* ; do for c in 30 40 50 ; do mkdir -p $f/c$c/illumina/bam && for r in $(seq -w 1 10) ; do minimap2 -t 4 -a -x sr ../reference/reference.fas $f/c$c/illumina/fastq/$f.c$c.illumina.r$r.fq.gz | samtools view -@ 4 -o $f/c$c/illumina/bam/$f.c$c.illumina.r$r.bam ; done ; done ; done
```

The individual Minimap2 command is as follows:

```bash
minimap2 -t THREADS -a -x sr REF_GENOME READS | samtools view -@ THREADS -o OUTPUT
```

* `-t THREADS` = Use `THREADS` threads
* `-a` = Output in the SAM format
* `-x sr` = Use the short-read preset
* `REF_GENOME` = FASTA file containing the reference genome
* `READS` = FASTQ file containing the reads
* `OUTPUT` = Output BAM file

# Oxford Nanopore Technologies (ONT)
We used [NanoSim-H v1.1.0.4]([https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm](https://github.com/karel-brinda/NanoSim-H/releases/tag/1.1.0.4)) to simulate ONT reads:

```bash
for f in lineage_* ; do for c in 10 30 50 ; do mkdir -p $f/c$c/ont/fasta && for r in $(seq -w 1 10) ; do nanosim-h -s $RANDOM -o $f/c$c/ont/fasta/$f.c$c.ont.r$r -n $(bc -l <<< "3.839 * $c" | numlist -ceil) $f/*.fas ; done ; done ; done
```

The individual NanoSim-H command is as follows:

```bash
nanosim-h -s RNG_SEED -o OUTPUT -n NUM_READS REF_GENOME
```

* `-s RNG_SEED` = Use `RNG_SEED` as the random number generation seed (`-s $RANDOM` picks a random value from the OS)
* `-o OUTPUT` = Output file prefix
* `-n NUM_READS` = Simulate `NUM_READS` reads
  * In its default settings, the average read length seems to be 7788.6321, so to obtain a coverage of roughly `C` for a genome of length `G`, one needs to simulate `G * C / 7788.6321` reads
  * For SARS-CoV-2, that's roughly 29903 * `C` / 7788.6321 = 3.839 * `C`
* `REF_GENOME` = FASTA file containing the reference genome

## Mapping
ONT reads were then mapped to the reference genome using Minimap2's ONT preset and piped to Samtools to convert to BAM:

```bash
for f in lineage_* ; do for c in 10 30 50 ; do mkdir -p $f/c$c/ont/bam && for r in $(seq -w 1 10) ; do minimap2 -t 4 -a -x map-ont ../reference/reference.fas $f/c$c/ont/fasta/$f.c$c.ont.r$r.fa.gz | samtools view -@ 4 -o $f/c$c/ont/bam/$f.c$c.ont.r$r.bam ; done ; done ; done
```

The individual Minimap2 command is as follows:

```bash
minimap2 -t THREADS -a -x map-ont REF_GENOME READS | samtools view -@ THREADS -o OUTPUT
```

* `-t THREADS` = Use `THREADS` threads
* `-a` = Output in the SAM format
* `-x map-ont` = Use the ONT preset
* `REF_GENOME` = FASTA file containing the reference genome
* `READS` = FASTQ file containing the reads
* `OUTPUT` = Output BAM file

# ViralConsensus
Consensus sequences were then called using ViralConsensus:

```bash
for f in */*/*/bam/*.bam ; do viral_consensus -r ../reference/reference.fas -i $f -o $(echo $f | sed 's/\.bam/.viralconsensus.fas/g' | sed 's/\/bam\//\/viralconsensus\//g') ; done
```

The individual ViralConsensus command is as follows:

```bash
viral_consensus -r REF_GENOME -i INPUT_BAM -o OUTPUT_FAS
```

* `-r REF_GENOME` = Use the reference genome in the FASTA file called `REF_GENOME`
* `-i INPUT_BAM` = Input from the BAM file called `INPUT_BAM`
* `-o OUTPUT_FAS` = Output to a FASTA file called `OUTPUT_FAS`

# iVar Pipeline
Consensus sequences were also called using the iVar pipeline.

## Sorting the BAMs
The Minimap2-mapped BAMs were sorted using Samtools:

```bash
for f in */*/*/bam/*.bam ; do samtools sort -@ 8 -o $(echo $f | sed 's/\.bam$/.sorted.bam/g') $f ; done
```

The individual Samtools command is as follows:

```bash
samtools view -@ THREADS -o OUTPUT_SORTED_BAM INPUT_BAM
```

* `-@ THREADS` = Use `THREADS` threads
* `-o OUTPUT_SORTED_BAM` = Output to a sorted BAM file called `OUTPUT_SORTED_BAM`
* `INPUT_BAM` = Input from the BAM file called `INPUT_BAM`

## Generating the Pile-ups
Pile-up files were calculated from the sorted BAMs using Samtools:

```bash
for f in */*/*/bam/*.sorted.bam ; do samtools mpileup -A -aa -d 0 -Q 0 --reference ../reference/reference.fas $f | pigz -9 -p 8 > $(echo $f | sed 's/\.bam$/.pileup.txt.gz/g' | sed 's/\/bam\//\/pileup\//g') ; done
```

The individual command is as follows:

```bash
samtools mpileup -A -aa -d 0 -Q 0 --reference REF_GENOME INPUT_SORTED_BAM | pigz -9 -p THREADS > OUTPUT_PILEUP
```

* `-A` = Count orphan reads
* `-aa` = Output absolutely all positions of the reference
* `-d 0` = No maximum depth
* `-Q 0` = Use a minimum base quality threshold of 0
* `--reference REF_GENOME` = Use the reference genome in the FASTA file called `REF_GENOME`
* `-9` = Use maximum GZIP compression
* `-p THREADS` = Compress using `THREADS` threads
* `OUTPUT_PILEUP_GZ` = The output gzip-compressed pile-up file

## Calling Consensus using iVar Consensus
Consensus sequences were called using iVar:

```bash
for f in */*/*/*/*.pileup.txt.gz ; do zcat $f | ivar consensus -p $(echo $f | sed 's/\/pileup\//\/ivarconsensus\//g' | sed 's/\.txt\.gz$/.ivar/g') -m 10 -n N -t 0.5 ; done
```

The individual command is as follows:

```bash
zcat INPUT_PILEUP_GZ | ivar consensus -p OUTPUT_PREFIX -m MIN_DEPTH -n AMBIG -t MIN_FREQ
```

* `INPUT_PILEUP_GZ` = The input gzip-compressed pile-up file
* `-p OUTPUT_PREFIX` = Use `OUTPUT_PREFIX` as the output file prefix
* `-m MIN_DEPTH` = Only call non-ambiguous bases in positions with a coverage of at least `MIN_DEPTH`
* `-n AMBIG` = Use `AMBIG` as the symbol representing ambiguous positions
* `-t MIN_FREQ` = Only call non-ambiguous bases in positions with a max base frequency of at least `MIN_FREQ`

# Pairwise Align to True Genome
Consensus sequences were pairwise-aligned to the true genome sequences using MAFFT:

```bash
for f in */*/*/*/*.viralconsensus.fas */*/*/*/*.ivar.fa ; do mafft <(cat $(echo $f | cut -d'/' -f1)/*.fas $f) | fasta1ln.py > $(echo $f | rev | cut -d'.' -f2- | rev).aln ; done
```

The individual command is as follows:

```bash
mafft <(cat TRUE_GENOME INPUT_CONSENSUS) | fasta1ln.py > OUTPUT_ALIGNMENT
```

* `TRUE_GENOME` = FASTA file containing the true genome
* `INPUT_CONSENSUS` = FASTA file containing the ViralConsensus/iVar consensus genome sequence
* `OUTPUT_ALIGNMENT` = FASTA file containing output pairwise alignment between `TRUE_GENOME` and `INPUT_CONSENSUS`
