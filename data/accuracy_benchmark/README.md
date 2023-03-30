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
for f in lineage_* ; do for c in 10 30 50 ; do mkdir -p $f/c$c/illumina/fastq && for r in $(seq -w 1 10) ; do art_illumina -rs $RANDOM -q -na -ss HS20 -l 100 -f $c -i $f/*.fas -o $f/c$c/illumina/fastq/$f.c$c.illumina.r$r ; done ; done ; done
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
TODO
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
for f in lineage_* ; do for c in 10 30 50 ; do mkdir -p $f/c$c/ont/fasta && for r in $(seq -w 1 10) ; do nanosim-h -s $RANDOM -o $f/c$c/ont/fasta/$f.c$c.illumina.r$r -n $(bc -l <<< "3.839 * $c") $f/*.fas ; done ; done ; done
```

The individual ART command is as follows:

```bash
nanosim-h -s RNG_SEED -o OUTPUT -n NUM_READS REF_GENOME
```

* `-s RNG_SEED` = Use `RNG_SEED` as the random number generation seed (`-s $RANDOM` picks a random value from the OS)
* `-o OUTPUT` = Output file prefix
* `-n NUM_READS` = Simulate `NUM_READS` reads
  * In its default settings, the average read length seems to be 7788.6321, so to obtain a coverage of roughly `C` for a genome of length `G`, one needs to simulate `G * C / 7788.6321` reads
  * For SARS-CoV-2, that's roughly 29903 * `C` / 7788.6321 = 3.839 * `C`
* `REF_GENOME` = FASTA file containing the reference genome
