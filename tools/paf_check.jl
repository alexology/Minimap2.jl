# Compare PAF records (no base-level alignment) with minimap2 (development check).
#   julia --project=. tools/paf_check.jl <minimap2> <preset> <reference.fa> <queries.fq>
using Minimap2
include(joinpath(@__DIR__, "..", "test", "fastx.jl"))

minimap2, preset_name, reference, queries = ARGS[1:4]
with_cigar = "-c" in ARGS[5:end]
refs = read_fastx(reference)
io, mo = preset(preset_name)
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
with_cigar && (mo.flag |= Minimap2.F_CIGAR | Minimap2.F_OUT_CG)
expected = split(replace(read(pipeline(`$minimap2 -x $preset_name -t 1 $(with_cigar ? ["-c"] : String[]) $reference $queries`; stderr=devnull), String), "\r" => ""), '\n'; keepempty=false)
got = String[]
elapsed = @elapsed for (name, seq) in read_fastx(queries)
    regions, rep_len = map_query(index, seq, mo; qname=name)
    for r in regions
        push!(got, paf_line(index, name, length(seq), r, rep_len; flag=mo.flag))
    end
end
diffs = count(i -> i > length(expected) || got[i] != expected[i], eachindex(got))
println("$preset_name $(with_cigar ? "-c" : ""): $(round(elapsed, digits=2)) s; $(length(got)) records vs $(length(expected)); differing records: $diffs")
if diffs > 0
    i = findfirst(i -> i > length(expected) || got[i] != expected[i], eachindex(got))
    println("  got:      ", got[i]); println("  expected: ", i <= length(expected) ? expected[i] : "-")
end
