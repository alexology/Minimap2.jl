using Minimap2
using Test
using CodecZlib

include("fastx.jl")
include("cases.jl")

reference_lines(name) = split(String(transcode(GzipDecompressor, read(joinpath(REFERENCE, name)))), '\n'; keepempty=false)

function setup(refs, preset_name)
    references = read_fastx(joinpath(DATA, refs))
    io, mo = preset(preset_name)
    index = MinimapIndex(last.(references), io; names=first.(references))
    Minimap2.update_options!(mo, index)
    return index, mo
end

@testset "primitives" begin
    @test Minimap2._mg_log2(8.0f0) ≈ 3.0 atol = 0.01
    v = [MM128(UInt64(x), UInt64(i)) for (i, x) in enumerate(rand(UInt64(0):UInt64(5), 500))]
    sorted = Minimap2.radix_sort_128x!(copy(v))
    @test issorted([a.x for a in sorted])
    @test sort([a.y for a in sorted]) == sort([a.y for a in v])
    minimizers = sketch!(MM128[], Vector{UInt8}(codeunits("ACGTACGTTAGCATCGATCGACTAGCTAGCATCGACTGACTAGCTACG")), 5, 11, 0, false)
    @test !isempty(minimizers)
    @test_throws ArgumentError preset("unknown")
    # minimap2's int seed wraps to 32 bits: -1 and 2^32 - 1 give the same tie-breaking hash.
    _, mo = preset("map-ont")
    mo.seed = -1; h = Minimap2._read_hash("read1", 1000, mo)
    mo.seed = 0xffffffff; @test Minimap2._read_hash("read1", 1000, mo) == h
end

@testset "index" begin
    # Above 1 Mb the index is built by several tasks; it must not depend on their number.
    rng_state = UInt64(12345)
    nextbase() = (rng_state = rng_state * 6364136223846793005 + 1442695040888963407; "ACGT"[(rng_state >> 62) + 1])
    unit = [String([nextbase() for _ in 1:2000]) for _ in 1:40]
    seqs = [unit[i % 40 + 1][1:(500 + 37i % 1500)] * unit[(7i) % 40 + 1] for i in 1:400]
    seqs[3] = ""; seqs[5] = "ACGTN" ^ 200
    for (k, w, hpc) in ((15, 10, false), (19, 19, true))
        indexes = [MinimapIndex(seqs; k, w, hpc, threads) for threads in (1, 4)]
        for f in (:keys, :starts, :positions, :directory, :directory_shift, :codes, :offsets, :lengths)
            @test getfield(indexes[1], f) == getfield(indexes[2], f)
        end
        table = Dict{UInt64,Vector{UInt64}}()
        for (rid, s) in enumerate(seqs)
            isempty(s) && continue
            for m in sketch!(MM128[], Vector{UInt8}(codeunits(s)), w, k, rid - 1, hpc)
                push!(get!(Vector{UInt64}, table, m.x >> 8), m.y)
            end
        end
        index = indexes[2]
        @test length(index.keys) == length(table)
        @test all(sort(v) == Minimap2.index_get(index, key) for (key, v) in table)
        @test isempty(Minimap2.index_get(index, typemax(UInt64)))
        @test all(isempty(Minimap2.index_get(index, key + 1)) for key in index.keys if !haskey(table, key + 1))
        counts = sort!([length(v) for v in values(table)])
        @test max_occurrence(index, 2.0f-4) == counts[unsafe_trunc(UInt32, (1.0 - Float64(2.0f-4)) * length(counts)) + 1] + 1
        # Any real fraction is accepted and rounded to Float32 like minimap2's float option.
        @test max_occurrence(index, 2.0e-4) == max_occurrence(index, 2.0f-4)
        @test max_occurrence(index, 0.01) == max_occurrence(index, 0.01f0)
        @test max_occurrence(index, 0) == typemax(Int32)
        @test max_occurrence(index, 1) == counts[1] + 1
        @test_throws ArgumentError max_occurrence(index, 1.5)
        @test_throws ArgumentError max_occurrence(index, NaN)
    end
    small = MinimapIndex(["ACGTACGTTAGCATCGATCGACTAGC"]; k=5, w=3)
    @test !isempty(small.keys) && Minimap2.reference_codes(small, 0, 0, 4) == UInt8[0, 1, 2, 3]
    @test_throws ArgumentError MinimapIndex(["ACGT"]; threads=0)
end

@testset "seeds and chains ($dataset)" for (dataset, refs, reads) in DATASETS
    index, mo = setup(refs, "map-ont")
    expected = reference_lines(debug_name(dataset, "map-ont"))
    got = String[]
    for (name, seq) in read_fastx(joinpath(DATA, reads))
        push!(got, "QR\t$name\t0\t$(length(seq))")
        anchors, rep_len = map_anchors(index, seq, mo; qname=name)
        push!(got, "RS\t$rep_len")
        for (i, a) in enumerate(anchors)
            diff = i == 1 ? 0 : (Int32(a.y & 0xffffffff) - Int32(anchors[i - 1].y & 0xffffffff)) - (Int32(a.x & 0xffffffff) - Int32(anchors[i - 1].x & 0xffffffff))
            push!(got, "SD\t$(index.names[Int((a.x << 1) >> 33) + 1])\t$(Int32(a.x & 0xffffffff))\t$("+-"[Int(a.x >> 63) + 1])\t$(Int32(a.y & 0xffffffff))\t$(Int((a.y >> 32) & 0xff))\t$diff")
        end
        regions, a, _, _ = map_chains(index, seq, mo; qname=name)
        for (j, r) in enumerate(regions), i in r.as:(r.as + r.cnt - 1)
            x = a[i + 1]
            diff = i == r.as ? 0 : (Int32(x.y & 0xffffffff) - Int32(a[i].y & 0xffffffff)) - (Int32(x.x & 0xffffffff) - Int32(a[i].x & 0xffffffff))
            push!(got, "CN\t$(j - 1)\t$(index.names[Int((x.x << 1) >> 33) + 1])\t$(Int32(x.x & 0xffffffff))\t$("+-"[Int(x.x >> 63) + 1])\t$(Int32(x.y & 0xffffffff))\t$(Int((x.y >> 32) & 0xff))\t$diff")
        end
    end
    # minimap2 prints chains after the seeds of the same read; compare per kind.
    for kind in ("QR\t", "RS\t", "SD\t", "CN\t")
        @test filter(l -> startswith(l, kind), got) == filter(l -> startswith(l, kind), expected)
    end
end

@testset "PAF $dataset $preset_name $(cigar ? "with CIGAR" : "")" for (dataset, refs, reads) in DATASETS,
        preset_name in PRESETS, cigar in (false, true)
    index, mo = setup(refs, preset_name)
    cigar && (mo.flag |= Minimap2.F_CIGAR | Minimap2.F_OUT_CG)
    got = String[]
    for (name, seq) in read_fastx(joinpath(DATA, reads))
        regions, rep_len = map_query(index, seq, mo; qname=name)
        append!(got, paf_line.(Ref(index), name, length(seq), regions, rep_len; flag=mo.flag))
    end
    @test got == reference_lines(paf_name(dataset, preset_name, cigar))
end

@testset "concurrent mapping ($dataset)" for (dataset, refs, reads) in DATASETS
    # map_query keeps no shared state: reads mapped by several tasks give the same records.
    index, mo = setup(refs, "map-ont")
    mo.flag |= Minimap2.F_CIGAR | Minimap2.F_OUT_CG
    queries = read_fastx(joinpath(DATA, reads))
    records = Vector{Vector{String}}(undef, length(queries))
    Threads.@threads for i in eachindex(queries)
        name, seq = queries[i]
        regions, rep_len = map_query(index, seq, mo; qname=name)
        records[i] = paf_line.(Ref(index), name, length(seq), regions, rep_len; flag=mo.flag)
    end
    @test reduce(vcat, records) == reference_lines(paf_name(dataset, "map-ont", true))
end
