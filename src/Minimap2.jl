"""
    Minimap2

Native Julia port of minimap2 2.30 for single-segment, unspliced long-read mapping.
It gives the same seeds, chains, mapping qualities, coordinates and CIGARs as
minimap2. Main entry points: `preset`, `MinimapIndex`,
`map_query` and `paf_line`.
"""
module Minimap2

using Printf: @sprintf

include("options.jl")
include("sort.jl")
include("sketch.jl")
include("index.jl")
include("seed.jl")
include("krmq.jl")
include("chain.jl")
include("regions.jl")
include("select.jl")
include("ksw2.jl")
include("ksw_extd2.jl")
include("ksw_ll.jl")
include("align.jl")
include("map.jl")

export IndexOptions, MapOptions, preset, MM128, sketch!, MinimapIndex, max_occurrence
export Region, AlignmentExtra, map_anchors, map_chains, map_query, paf_line

end
