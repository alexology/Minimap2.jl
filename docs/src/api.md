# API

```@meta
CurrentModule = Minimap2
```

```@docs
Minimap2
```

## Options

```@docs
IndexOptions
MapOptions
preset
update_options!
```

Mapping flags are the `MapOptions.flag` bits of minimap2 (`MM_F_*` in `minimap.h`),
available as unexported constants such as `Minimap2.F_CIGAR` (base-level alignment)
and `Minimap2.F_OUT_CG` (`cg:Z` tag).

## Index

```@docs
MM128
sketch!
MinimapIndex
max_occurrence
```

## Mapping

```@docs
Region
AlignmentExtra
map_anchors
map_chains
map_query
paf_line
```
