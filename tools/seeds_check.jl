# Compare anchors with `minimap2 --print-seeds` (development check).
#   julia --project=. tools/seeds_check.jl <minimap2> <preset> <reference.fa> <queries.fq>
using Minimap2

function read_fastx(path)
    records = Tuple{String,Vector{UInt8}}[]
    lines = readlines(path)
    i = 1
    while i <= length(lines)
        if startswith(lines[i], '>')
            name = split(lines[i][2:end])[1]; seq = UInt8[]; i += 1
            while i <= length(lines) && !startswith(lines[i], '>')
                append!(seq, codeunits(lines[i])); i += 1
            end
            push!(records, (name, seq))
        else
            push!(records, (split(lines[i][2:end])[1], Vector{UInt8}(codeunits(lines[i + 1])))); i += 4
        end
    end
    return records
end

minimap2, preset_name, reference, queries = ARGS
refs = read_fastx(reference)
io, mo = preset(preset_name)
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
err = IOBuffer()
run(pipeline(`$minimap2 -x $preset_name -t 1 --print-seeds $reference $queries`; stdout=devnull, stderr=err))
blocks = Vector{Vector{String}}(); current = String[]
for line in split(replace(String(take!(err)), "\r" => ""), '\n')
    if startswith(line, "RS\t")
        global current = String[line]; push!(blocks, current)
    elseif startswith(line, "SD\t")
        push!(current, line)
    end
end
bad = 0
for (qi, (name, seq)) in enumerate(read_fastx(queries))
    anchors, rep_len = map_anchors(index, seq, mo; qname=name)
    got = ["RS\t$rep_len"]
    for (i, a) in enumerate(anchors)
        diff = i == 1 ? 0 : (Int32(a.y & 0xffffffff) - Int32(anchors[i - 1].y & 0xffffffff)) - (Int32(a.x & 0xffffffff) - Int32(anchors[i - 1].x & 0xffffffff))
        push!(got, "SD\t$(index.names[Int((a.x << 1) >> 33) + 1])\t$(Int32(a.x & 0xffffffff))\t$("+-"[Int(a.x >> 63) + 1])\t$(Int32(a.y & 0xffffffff))\t$(Int((a.y >> 32) & 0xff))\t$diff")
    end
    if got != blocks[qi]
        global bad += 1
        bad <= 2 && println("query $qi: $(length(got)) vs $(length(blocks[qi])) lines; first difference at ", findfirst(i -> i > length(got) || i > length(blocks[qi]) || got[i] != blocks[qi][i], 1:max(length(got), length(blocks[qi]))))
    end
end
println("mid_occ = $(mo.mid_occ); queries with different anchors: $bad / $(length(blocks))")
