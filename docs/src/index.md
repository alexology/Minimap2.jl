# Minimap2.jl

Native Julia port of [minimap2](https://github.com/lh3/minimap2) 2.30, written to
give the same results as minimap2: same seeds, chains, mapping qualities,
coordinates and CIGARs. It was started for
[Savont.jl](https://github.com/alexology/Savont.jl), whose upstream (savont 0.8.0)
calls minimap2 2.30 through the Rust crate, and is meant to be shared with other
packages of this suite.

This package is part of an experimental effort to use AI to port bioinformatics
and quantitative ecology packages to Julia, with the goal of creating a unified
environment for high performance scientific computing.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/Minimap2.jl")
```

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

Ties between equal chains are broken by a hash seeded with the read name, as in
minimap2: pass `qname` to [`map_query`](@ref) to reproduce the command line, or
leave it out to match library calls without names (minimap2's `--no-hash-name`,
the Rust crate called with `None`). `MapOptions.seed` also enters the hash and,
like minimap2's `int` seed, wraps to 32 bits.

[`map_query`](@ref) keeps no shared state, so different reads can be mapped from
different threads with the same index and options.

## Status

Implemented for single-segment, unspliced mapping:

- presets `map-ont` (`lr`), `lr:hq` (`map-hifi`, `map-ccs`), `map-pb`
  (homopolymer-compressed seeds), `asm5/10/20` and `sr` options;
- minimizer sketches, in-memory index (minimizers sorted by hash in flat arrays,
  built by several tasks with `threads`; the result does not depend on their
  number), seed collection with minimap2's high-occurrence selection;
- DP and RMQ chaining (minimap2's AVL range-minimum tree is ported line by line,
  since ties depend on the tree shape), long-join re-chaining;
- primary/secondary selection, divergence estimate, mapping quality;
- base-level alignment: `ksw_extd2` (dual affine gaps) and the striped
  `ksw_ll_i16`, emulated cell by cell with the arithmetic of their SSE code, gap
  filling, Z-drop, inversion rescue, CIGAR left-alignment;
- PAF records ([`paf_line`](@ref)), including `cg:Z`.

Not implemented: single affine gap alignment (`ksw_extz2`, used only when the two
gap penalties coincide), spliced alignment, the short-read path (`sr` mapping),
paired ends, SAM output, SDUST masking, saving the index (`.mmi`) and
multi-threaded mapping of a batch of reads.

minimap2's in-place radix sort is not stable; it is ported literally so that equal
keys come out in the same order.

## License and attribution

minimap2 is by Heng Li and contributors, MIT license; the package `LICENSE`
reproduces it.

- Li H. *Minimap2: pairwise alignment for nucleotide sequences.* Bioinformatics
  34:3094–3100 (2018).
- Li H. *New strategies to improve minimap2 alignment accuracy.* Bioinformatics
  37:4572–4574 (2021).
