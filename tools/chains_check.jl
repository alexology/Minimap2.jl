# Compare chains with `minimap2 --print-chains` (development check).
#   julia --project=. tools/chains_check.jl <minimap2> <preset> <reference.fa> <queries.fq> [extra minimap2 options]
using Minimap2
include(joinpath(@__DIR__, "..", "test", "fastx.jl"))

minimap2, preset_name, reference, queries = ARGS[1:4]
extra = ARGS[5:end]
refs = read_fastx(reference)
io, mo = preset(preset_name)
"--rmq" in extra && (mo.flag |= Minimap2.F_RMQ)
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
err = IOBuffer()
run(pipeline(`$minimap2 -x $preset_name -t 1 --print-qname --print-chains $extra $reference $queries`; stdout=devnull, stderr=err))
# Chains of each query follow its QR line.
blocks = Vector{Vector{String}}(); current = String[]
for line in split(replace(String(take!(err)), "\r" => ""), '\n')
    if startswith(line, "QR\t")
        global current = String[]; push!(blocks, current)
    elseif startswith(line, "CN\t")
        push!(current, line)
    end
end
bad = 0; n_chains = 0
for (qi, (name, seq)) in enumerate(read_fastx(queries))
    regions, a, _, _ = map_chains(index, seq, mo; qname=name)
    got = String[]
    for (j, r) in enumerate(regions), i in r.as:(r.as + r.cnt - 1)
        x = a[i + 1]
        diff = i == r.as ? 0 : (Int32(x.y & 0xffffffff) - Int32(a[i].y & 0xffffffff)) - (Int32(x.x & 0xffffffff) - Int32(a[i].x & 0xffffffff))
        push!(got, "CN\t$(j - 1)\t$(index.names[Int((x.x << 1) >> 33) + 1])\t$(Int32(x.x & 0xffffffff))\t$("+-"[Int(x.x >> 63) + 1])\t$(Int32(x.y & 0xffffffff))\t$(Int((x.y >> 32) & 0xff))\t$diff")
    end
    global n_chains += length(regions)
    if got != blocks[qi]
        global bad += 1
        bad <= 2 && println("query $qi ($name): $(length(got)) vs $(length(blocks[qi])) lines; first difference at ", findfirst(i -> i > length(blocks[qi]) || got[i] != blocks[qi][i], 1:length(got)))
    end
end
println("$preset_name $(join(extra, ' ')): queries with different chains: $bad / $(length(blocks)) ($n_chains chains)")
