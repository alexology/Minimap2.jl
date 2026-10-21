# Profile mapping with base-level alignment (development tool).
#   julia --project=. tools/profile_map.jl <reference.fa> <queries.fq> [n_queries]
using Minimap2, Profile
include(joinpath(@__DIR__, "..", "test", "fastx.jl"))
refs = read_fastx(ARGS[1]); queries = read_fastx(ARGS[2])
n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : length(queries)
io, mo = preset("map-ont")
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
mo.flag |= Minimap2.F_CIGAR
run_all() = for (name, seq) in queries[1:n]; map_query(index, seq, mo; qname=name); end
map_query(index, queries[1][2], mo; qname=queries[1][1])
t = @elapsed run_all()
alloc = @allocated run_all()
println("$(n) queries: $(round(t, digits=2)) s, $(round(alloc / 2^30, digits=2)) GiB allocated")
Profile.clear(); @profile run_all()
Profile.print(IOContext(stdout, :displaysize => (200, 200)); format=:flat, sortedby=:count, mincount=20, noisefloor=2)
