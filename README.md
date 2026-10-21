# Minimap2.jl

Native Julia port of [minimap2](https://github.com/lh3/minimap2) 2.30 (`UPSTREAM.toml`),
written to give the same results as minimap2: same seeds, chains, mapping
qualities, coordinates and CIGARs. It was started for
[Savont.jl](../Savont.jl), whose upstream (savont 0.8.0) calls minimap2 2.30
through the Rust crate, and is meant to be shared with other packages of this
suite.

## Status

Implemented for single-segment, unspliced mapping:

- presets `map-ont` (`lr`), `lr:hq` (`map-hifi`, `map-ccs`), `map-pb` (homopolymer-compressed seeds), `asm5/10/20` and `sr` options;
- minimizer sketches, in-memory index, seed collection with minimap2's high-occurrence selection;
- DP and RMQ chaining (minimap2's AVL range-minimum tree is ported line by line, since ties depend on the tree shape), long-join re-chaining;
- primary/secondary selection, divergence estimate, mapping quality;
- base-level alignment: `ksw_extd2` (dual affine gaps) and the striped `ksw_ll_i16`, emulated cell by cell with the arithmetic of their SSE code, gap filling, Z-drop, inversion rescue, CIGAR left-alignment;
- PAF records (`paf_line`), including `cg:Z`.

Not implemented: single affine gap alignment (`ksw_extz2`, used only when the two gap penalties coincide), spliced alignment, the short-read path (`sr` mapping), paired ends, SAM output, SDUST masking and multi-threading (map different reads on different threads).

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

`test/runtests.jl` compares, on three data sets (amplicons, reads with 2% substitutions, and reads with inversions, chimeras, long indels, runs of `N` and 15% errors), the anchors and chains printed by `minimap2 --print-seeds --print-chains` (`map-ont`) and the PAF records with and without CIGAR for `map-ont`, `lr:hq` and `map-pb`. All records are identical. The expected outputs are stored in `test/reference` and regenerated with `julia tools/make_reference.jl`, which runs minimap2 built from the pinned checkout:

```sh
git clone --branch v2.30 https://github.com/lh3/minimap2 reference/minimap2
make -C reference/minimap2        # on Windows: gmake CC=gcc (MinGW)
```

`tools/` also holds development checks (`seeds_check.jl`, `chains_check.jl`, `paf_check.jl`) that run minimap2 on any FASTA/FASTQ input.

## Performance

Single thread, 3,190 reads of 1.2 kb against 50 similar references: 25 s with base-level alignment (minimap2: 15 s) and 16 s without (minimap2: 12 s). The alignment kernel runs at 0.4–1.3 ns per DP cell. `Pkg.test()` runs with bounds checking forced on and is slower.

## License and attribution

minimap2 is by Heng Li and contributors, MIT license; `LICENSE` reproduces it. Li H. *Minimap2: pairwise alignment for nucleotide sequences.* Bioinformatics 34:3094–3100 (2018); Li H. *New strategies to improve minimap2 alignment accuracy.* Bioinformatics 37:4572–4574 (2021).
