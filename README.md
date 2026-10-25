# Minimap2.jl

Native Julia port of [minimap2](https://github.com/lh3/minimap2) 2.30 (`UPSTREAM.toml`),
written to give the same results as minimap2: same seeds, chains, mapping
qualities, coordinates and CIGARs. It was started for
[Savont.jl](https://github.com/alexology/Savont.jl), whose upstream (savont 0.8.0) calls minimap2 2.30
through the Rust crate, and is meant to be shared with other packages of this
suite.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

## Install from GitHub

This package is not registered in the Julia General registry. Install it
from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/Minimap2.jl")
```

## Status

Implemented for single-segment, unspliced mapping:

- presets `map-ont` (`lr`), `lr:hq` (`map-hifi`, `map-ccs`), `map-pb` (homopolymer-compressed seeds), `asm5/10/20` and `sr` options;
- minimizer sketches, in-memory index (minimizers sorted by hash in flat arrays, built by several tasks with `threads`; the result does not depend on their number), seed collection with minimap2's high-occurrence selection;
- DP and RMQ chaining (minimap2's AVL range-minimum tree is ported line by line, since ties depend on the tree shape), long-join re-chaining;
- primary/secondary selection, divergence estimate, mapping quality;
- base-level alignment: `ksw_extd2` (dual affine gaps) and the striped `ksw_ll_i16`, emulated cell by cell with the arithmetic of their SSE code, gap filling, Z-drop, inversion rescue, CIGAR left-alignment;
- PAF records (`paf_line`), including `cg:Z`.

Not implemented: single affine gap alignment (`ksw_extz2`, used only when the two gap penalties coincide), spliced alignment, the short-read path (`sr` mapping), paired ends, SAM output, SDUST masking, saving the index (`.mmi`) and multi-threaded mapping (`map_query` can be called on different reads from different threads).

minimap2's in-place radix sort is not stable; it is ported literally so that equal keys come out in the same order.

## Usage

```julia
using Minimap2
index_options, map_options = preset("map-ont")
index = MinimapIndex(reference_sequences, index_options; names=reference_names)
Minimap2.update_options!(map_options, index)
map_options.flag |= Minimap2.F_CIGAR              # base-level alignment
regions, rep_len = map_query(index, read_sequence, map_options; qname="read1")
for r in regions
    println(paf_line(index, "read1", length(read_sequence), r, rep_len; flag=map_options.flag))
end
```

## Validation

`test/runtests.jl` compares, on four data sets (amplicons, reads with 2% substitutions, reads with inversions, chimeras, long indels, runs of `N` and 15% errors, and a read whose inverted part ends in the padding of the striped query profile), the anchors and chains printed by `minimap2 --print-seeds --print-chains` (`map-ont`) and the PAF records with and without CIGAR for `map-ont`, `lr:hq` and `map-pb`. All records are identical. The expected outputs are stored in `test/reference` and regenerated with `julia tools/make_reference.jl`, which runs minimap2 built from the pinned checkout:

```sh
git clone --branch v2.30 https://github.com/lh3/minimap2 reference/minimap2
make -C reference/minimap2        # on Windows: gmake CC=gcc (MinGW)
```

`tools/` also holds development checks (`seeds_check.jl`, `chains_check.jl`, `paf_check.jl`) that run minimap2 on any FASTA/FASTQ input.

On 20,000 real ONT reads (~2.8 kb) against 32 amplicons, `paf_check.jl` gives the same 62,292 (`map-ont`) and 54,879 (`lr:hq`) records without base-level alignment, with and without `--no-hash-name`. With base-level alignment the only differences (8 of 62,092 and 6 of 54,374 records) are `de:f:0.0312` against `de:f:0.0313` for a divergence of exactly 0.03125: the MinGW C runtime used by the Windows build of minimap2 rounds the tie up, while glibc and Julia round it to even. On 96 synthetic reads with inversions, deletions, insertions and chimeras (1% substitutions), the records are identical for `map-ont`, `lr:hq`, `map-pb` and `asm5/10/20`, with and without base-level alignment.

Ties between equal chains are broken by a hash seeded with the read name, as in minimap2: pass `qname` to `map_query` to reproduce the command line, or leave it out to match library calls without names (minimap2's `--no-hash-name`, the Rust crate called with `None`). `MapOptions.seed` also enters the hash and, like minimap2's `int` seed, wraps to 32 bits.

## Performance

Index of 50,000 references of 1,450 random bases (72.5 Mb): 2.7 s on one thread, 1.1 s on four, 0.42 GiB.

Single thread, 3,190 reads of 1.2 kb against 50 similar references: 25 s with base-level alignment (minimap2: 15 s) and 16 s without (minimap2: 12 s). The alignment kernel runs at 0.4–1.3 ns per DP cell. `Pkg.test()` runs with bounds checking forced on and is slower.

## License and attribution

minimap2 is by Heng Li and contributors, MIT license; `LICENSE` reproduces it. Li H. *Minimap2: pairwise alignment for nucleotide sequences.* Bioinformatics 34:3094–3100 (2018); Li H. *New strategies to improve minimap2 alignment accuracy.* Bioinformatics 37:4572–4574 (2021).
