# Compare PAF records with minimap2 (development check).
#   julia --project=. tools/paf_check.jl <minimap2> <preset> <reference.fa> <queries.fq> [-c] [--no-hash-name] [--seed N]
# `-c` adds base-level alignment; `--no-hash-name` seeds tie-breaking without the
# read name, as when minimap2 is called through its library without names;
# `--seed` changes the seed of that hash.
using Minimap2
include(joinpath(@__DIR__, "..", "test", "fastx.jl"))

minimap2, preset_name, reference, queries = ARGS[1:4]
options = ARGS[5:end]
with_cigar = "-c" in options
no_hash_name = "--no-hash-name" in options
seed_at = findfirst(==("--seed"), options)
refs = read_fastx(reference)
io, mo = preset(preset_name)
index = MinimapIndex(last.(refs), io; names=first.(refs))
Minimap2.update_options!(mo, index)
with_cigar && (mo.flag |= Minimap2.F_CIGAR | Minimap2.F_OUT_CG)
no_hash_name && (mo.flag |= Minimap2.F_NO_HASH_NAME)
seed_at === nothing || (mo.seed = parse(Int, options[seed_at + 1]))
extra = [with_cigar ? ["-c"] : String[]; no_hash_name ? ["--no-hash-name"] : String[];
         seed_at === nothing ? String[] : ["--seed", options[seed_at + 1]]]
expected = split(replace(read(pipeline(`$minimap2 -x $preset_name -t 1 $extra $reference $queries`; stderr=devnull), String), "\r" => ""), '\n'; keepempty=false)
got = String[]
elapsed = @elapsed for (name, seq) in read_fastx(queries)
    regions, rep_len = map_query(index, seq, mo; qname=name)
    for r in regions
        push!(got, paf_line(index, name, length(seq), r, rep_len; flag=mo.flag))
    end
end
# Records missing on either side count as differences.
differs(i) = i > length(got) || i > length(expected) || got[i] != expected[i]
records = 1:max(length(got), length(expected))
diffs = count(differs, records)
println("$preset_name $(join(extra, " ")): $(round(elapsed, digits=2)) s; $(length(got)) records vs $(length(expected)); differing records: $diffs")
if diffs > 0
    i = findfirst(differs, records)
    println("  got:      ", i <= length(got) ? got[i] : "-"); println("  expected: ", i <= length(expected) ? expected[i] : "-")
end
