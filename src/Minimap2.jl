module Minimap2

using Printf

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
