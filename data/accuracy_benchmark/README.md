To benchmark the accuracy of ViralConsensus, simulated reads from known genomes and then compared the ViralConsensus sequence vs. the true genome sequence. We selected a high-confidence representative genome from each of the following SARS-CoV-2 lineages:
* **Lineage B:** [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)
* **Lineage B.1.1.7:** [LC650844.1](https://www.ncbi.nlm.nih.gov/nuccore/LC650844.1)
* **Lineage BA.1:** [OQ523614.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ523614.1)
* **Lineage BA.2:** [OQ194009.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ194009.1)
* **Lineage XBB.1:** [OQ346068.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ346068.1)

# Illumina
We used [ART version MountRainier-2016-06-05](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to simulate Illumina reads:

```bash
for f in lineage_* ; do for c in 10 30 50 ; do mkdir -p $f/c$c/illumina/fastq && for r in $(seq -w 1 10) ; do art_illumina -rs $RANDOM -q -na -ss HS20 -l 100 -f $c -i $f/*.fas -o $f/c$c/illumina/$f.c$c.illumina.r$r ; done ; done ; done
```

The individual ART command is as follows:

```bash
art_illumina -rs RNG_SEED -q -na -ss HS20 -l 100 -f COVERAGE -i REF_GENOME -o OUTPUT
```

* `-rs RNG_SEED` = Use `RNG_SEED` as the random number generation seed (`-r $RANDOM` picks a random value from the OS)
* `-q` = Quiet mode
* `-na` = Don't output the alignment file (just the FASTQ)
* `-ss HS20` = Use the HiSeq 2000 error profile
* `-l 100` = Use a length of 100 bases for each read
* `-f COVERAGE` = Simulate at `COVERAGE` fold coverage (e.g. `-f 10` = 10X coverage)
* `-i REF_GENOME` = FASTA file containing the reference genome
* `-o OUTPUT` = Output file prefix
