# Indexing and mapping options with minimap2's presets (options.c).

# Mapping flags (minimap.h, MM_F_*).
const F_NO_DIAG = Int64(0x001)
const F_NO_DUAL = Int64(0x002)
const F_CIGAR = Int64(0x004)
const F_OUT_CG = Int64(0x020)
const F_SPLICE = Int64(0x080)
const F_SPLICE_FOR = Int64(0x100)
const F_SPLICE_REV = Int64(0x200)
const F_NO_LJOIN = Int64(0x400)
const F_SR = Int64(0x1000)
const F_FRAG_MODE = Int64(0x2000)
const F_NO_PRINT_2ND = Int64(0x4000)
const F_2_IO_THREADS = Int64(0x8000)
const F_SPLICE_FLANK = Int64(0x40000)
const F_FOR_ONLY = Int64(0x100000)
const F_REV_ONLY = Int64(0x200000)
const F_HEAP_SORT = Int64(0x400000)
const F_ALL_CHAINS = Int64(0x800000)
const F_NO_END_FLT = Int64(0x10000000)
const F_HARD_MLEVEL = Int64(0x20000000)
const F_RMQ = Int64(0x80000000)
const F_QSTRAND = Int64(0x100000000)
const F_NO_INV = Int64(0x200000000)
const F_NO_HASH_NAME = Int64(0x400000000)
const F_WEAK_PAIRING = Int64(0x4000000000)
const F_SR_RNA = Int64(0x8000000000)

"""
    IndexOptions(; k=15, w=10, hpc=false, bucket_bits=14)

Minimizer parameters of the index (minimap2 `mm_idxopt_t`): k-mer size `k`,
window `w` and homopolymer compression `hpc`.
"""
Base.@kwdef mutable struct IndexOptions
    k::Int = 15
    w::Int = 10
    hpc::Bool = false
    bucket_bits::Int = 14
end

"""
    MapOptions

Mapping options (minimap2 `mm_mapopt_t`), with minimap2's defaults. Field names
follow the C structure; see [`preset`](@ref) for the presets.
"""
Base.@kwdef mutable struct MapOptions
    flag::Int64 = 0
    seed::Int = 11
    sdust_thres::Int = 0
    max_qlen::Int = 0
    bw::Int = 500
    bw_long::Int = 20000
    max_gap::Int = 5000
    max_gap_ref::Int = -1
    max_frag_len::Int = 0
    max_chain_skip::Int = 25
    max_chain_iter::Int = 5000
    min_cnt::Int = 3
    min_chain_score::Int = 40
    chain_gap_scale::Float32 = 0.8f0
    chain_skip_scale::Float32 = 0.0f0
    rmq_size_cap::Int = 100000
    rmq_inner_dist::Int = 1000
    rmq_rescue_size::Int = 1000
    rmq_rescue_ratio::Float32 = 0.1f0
    mask_level::Float32 = 0.5f0
    mask_len::Int = typemax(Int32)
    pri_ratio::Float32 = 0.8f0
    best_n::Int = 5
    alt_drop::Float32 = 0.15f0
    a::Int = 2
    b::Int = 4
    q::Int = 4
    e::Int = 2
    q2::Int = 24
    e2::Int = 1
    transition::Int = 0
    sc_ambi::Int = 1
    noncan::Int = 0
    junc_bonus::Int = 0
    junc_pen::Int = 0
    zdrop::Int = 400
    zdrop_inv::Int = 200
    end_bonus::Int = -1
    min_dp_max::Int = 80
    min_ksw_len::Int = 200
    anchor_ext_len::Int = 20
    anchor_ext_shift::Int = 6
    max_clip_ratio::Float32 = 1.0f0
    rank_min_len::Int = 500
    rank_frac::Float32 = 0.9f0
    pe_ori::Int = 0
    pe_bonus::Int = 33
    jump_min_match::Int = 3
    mid_occ_frac::Float32 = 2.0f-4
    q_occ_frac::Float32 = 0.01f0
    min_mid_occ::Int = 10
    max_mid_occ::Int = 1000000
    mid_occ::Int = 0
    max_occ::Int = 0
    max_max_occ::Int = 4095
    occ_dist::Int = 500
    max_sw_mat::Int64 = 100000000
end

"""
    preset(name) -> (IndexOptions, MapOptions)

minimap2's presets (`mm_set_opt`, options `-x`): `"map-ont"` (also `"lr"`, the
default), `"lr:hq"` (`"map-hifi"`, `"map-ccs"`), `"map-pb"`, `"asm5"`,
`"asm10"`, `"asm20"` and `"sr"`.
"""
function preset(name::AbstractString)
    io = IndexOptions(); mo = MapOptions()
    if name == "lr" || name == "map-ont"
        # The defaults.
    elseif name == "map10k" || name == "map-pb"
        io.hpc = true; io.k = 19
    elseif name == "lr:hq" || name == "map-hifi" || name == "map-ccs"
        io.hpc = false; io.k = 19; io.w = 19
        mo.max_gap = 10000
        mo.min_mid_occ = 50; mo.max_mid_occ = 500
        if name != "lr:hq"
            mo.a = 1; mo.b = 4; mo.q = 6; mo.q2 = 26; mo.e = 2; mo.e2 = 1
            mo.min_dp_max = 200
        end
    elseif startswith(name, "asm")
        io.hpc = false; io.k = 19; io.w = 19
        mo.bw = 1000; mo.bw_long = 100000
        mo.max_gap = 10000
        mo.flag |= F_RMQ
        mo.min_mid_occ = 50; mo.max_mid_occ = 500
        mo.min_dp_max = 200
        mo.best_n = 50
        if name == "asm5"
            mo.a = 1; mo.b = 19; mo.q = 39; mo.q2 = 81; mo.e = 3; mo.e2 = 1; mo.zdrop = mo.zdrop_inv = 200
        elseif name == "asm10"
            mo.a = 1; mo.b = 9; mo.q = 16; mo.q2 = 41; mo.e = 2; mo.e2 = 1; mo.zdrop = mo.zdrop_inv = 200
        elseif name == "asm20"
            mo.a = 1; mo.b = 4; mo.q = 6; mo.q2 = 26; mo.e = 2; mo.e2 = 1; mo.zdrop = mo.zdrop_inv = 200
            io.w = 10
        else
            throw(ArgumentError("unknown preset: $name"))
        end
    elseif name == "short" || name == "sr"
        io.hpc = false; io.k = 21; io.w = 11
        mo.flag |= F_SR | F_FRAG_MODE | F_NO_PRINT_2ND | F_2_IO_THREADS | F_HEAP_SORT
        mo.pe_ori = 0 << 1 | 1
        mo.a = 2; mo.b = 8; mo.q = 12; mo.e = 2; mo.q2 = 24; mo.e2 = 1
        mo.zdrop = mo.zdrop_inv = 100
        mo.end_bonus = 10
        mo.max_frag_len = 800
        mo.max_gap = 100
        mo.bw = mo.bw_long = 100
        mo.pri_ratio = 0.5f0
        mo.min_cnt = 2
        mo.min_chain_score = 25
        mo.min_dp_max = 40
        mo.best_n = 20
        mo.mid_occ = 1000
        mo.max_occ = 5000
    else
        throw(ArgumentError("unknown or unsupported preset: $name"))
    end
    return io, mo
end

# mm_mapopt_update: derive the occurrence threshold from the index.
function update_options!(opt::MapOptions, index)
    if (opt.flag & F_SPLICE_FOR) != 0 || (opt.flag & F_SPLICE_REV) != 0
        opt.flag |= F_SPLICE
    end
    if opt.mid_occ <= 0
        opt.mid_occ = max_occurrence(index, opt.mid_occ_frac)
        opt.mid_occ < opt.min_mid_occ && (opt.mid_occ = opt.min_mid_occ)
        if opt.max_mid_occ > opt.min_mid_occ && opt.mid_occ > opt.max_mid_occ
            opt.mid_occ = opt.max_mid_occ
        end
    end
    opt.bw_long < opt.bw && (opt.bw_long = opt.bw)
    return opt
end
