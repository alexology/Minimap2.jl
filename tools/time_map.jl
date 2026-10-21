using Minimap2
include(joinpath(@__DIR__, "..", "test", "fastx.jl"))
refs = read_fastx(ARGS[1]); queries = read_fastx(ARGS[2])
io, mo = preset("map-ont")
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
"nocigar" in ARGS || (mo.flag |= Minimap2.F_CIGAR)
function run(index, mo, qs)
    times = Float64[]
    for (name, seq) in qs
        push!(times, @elapsed map_query(index, seq, mo; qname=name))
    end
    times
end
run(index, mo, queries[1:3])
t = run(index, mo, queries)
println("total ", round(sum(t), digits=2), " s; slowest ", round.(sort(t; rev=true)[1:5], digits=2), " median ", round(sort(t)[length(t)÷2], digits=3))
