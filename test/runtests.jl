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
