Run simulations:

```bash
for n in 100 1000 10000 100000 1000000 ; do mkdir -p n$n && for r in $(seq -w 1 10) ; do ../scripts/run_subsample.sh reads.fastq.gz $n reference.fas reference.fas.mmi primers.bed n$n/n$n.r$r.zip $RANDOM ; done ; done
```
