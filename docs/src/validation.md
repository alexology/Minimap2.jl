# Validation and performance

## Comparison with minimap2

The test suite compares, on four data sets (amplicons, reads with 2%
substitutions, reads with inversions, chimeras, long indels, runs of `N` and 15%
errors, and a read whose inverted part ends in the padding of the striped query
profile), the anchors and chains printed by
`minimap2 --print-seeds --print-chains` (`map-ont`) and the PAF records with and
without CIGAR for `map-ont`, `lr:hq` and `map-pb`. All records are identical.

The expected outputs are stored in `test/reference` and regenerated with
`julia tools/make_reference.jl`, which runs minimap2 built from the pinned
checkout (`UPSTREAM.toml`):

```sh
git clone --branch v2.30 https://github.com/lh3/minimap2 reference/minimap2
make -C reference/minimap2        # on Windows: gmake CC=gcc (MinGW)
```

`tools/` also holds development checks (`seeds_check.jl`, `chains_check.jl`,
`paf_check.jl`) that run minimap2 on any FASTA/FASTQ input.

On 20,000 real ONT reads (~2.8 kb) against 32 amplicons, `paf_check.jl` gives the
same 62,292 (`map-ont`) and 54,879 (`lr:hq`) records without base-level
alignment, with and without `--no-hash-name`. With base-level alignment the only
differences (8 of 62,092 and 6 of 54,374 records) are `de:f:0.0312` against
`de:f:0.0313` for a divergence of exactly 0.03125: the MinGW C runtime used by
the Windows build of minimap2 rounds the tie up, while glibc and Julia round it
to even. On 96 synthetic reads with inversions, deletions, insertions and
chimeras (1% substitutions), the records are identical for `map-ont`, `lr:hq`,
`map-pb` and `asm5/10/20`, with and without base-level alignment.

## Performance

Index of 50,000 references of 1,450 random bases (72.5 Mb): 2.7 s on one thread,
1.1 s on four, 0.42 GiB.

Single thread, 3,190 reads of 1.2 kb against 50 similar references: 25 s with
base-level alignment (minimap2: 15 s) and 16 s without (minimap2: 12 s). The
alignment kernel runs at 0.4–1.3 ns per DP cell. `Pkg.test()` runs with bounds
checking forced on and is slower.
