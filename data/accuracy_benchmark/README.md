To benchmark the accuracy of ViralConsensus, simulated reads from known genomes and then compared the ViralConsensus sequence vs. the true genome sequence. We selected a high-confidence representative genome from each of the following SARS-CoV-2 lineages:
* **Lineage B:** [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)
* **Lineage B.1.1.7:** [LC650844.1](https://www.ncbi.nlm.nih.gov/nuccore/LC650844.1)
* **Lineage BA.1:** [OQ523614.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ523614.1)
* **Lineage BA.2:** [OQ194009.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ194009.1)
* **Lineage XBB.1:** [OQ346068.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ346068.1)

# Illumina
We used [ART version MountRainier-2016-06-05](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to simulate Illumina reads:

```bash
art_illumina -ss HS25 -l 100 -f COVERAGE -i REF_GENOME -o OUTPUT
```

* `-ss HS20` = Use the HiSeq 2000 error profile
* `-l 100` = Use a length of 100 bases for each read
* `-f COVERAGE` = Simulate at `COVERAGE` fold coverage (e.g. `-f 10` = 10X coverage)
* `-i REF_GENOME` = FASTA file containing the reference genome
* `-o OUTPUT` = Output file
