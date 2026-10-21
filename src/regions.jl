# Mapping regions (mm_reg1_t) and their construction from chains (hit.c).

const PARENT_UNSET = -1
const PARENT_TMP_PRI = -2

"""
    AlignmentExtra

Base-level alignment of a region (minimap2 `mm_extra_t`): DP score, maximum
scores, number of ambiguous bases and the CIGAR as `length << 4 | op` words with
ops `MIDNSHP=XB` numbered 0–9.
"""
mutable struct AlignmentExtra
    dp_score::Int32
    dp_max::Int32
    dp_max2::Int32
    dp_max0::Int32
    n_ambi::UInt32
    trans_strand::UInt32
    cigar::Vector{UInt32}
end

"""
    Region

A mapping (minimap2 `mm_reg1_t`), with 0-based half-open query (`qs`, `qe`) and
reference (`rs`, `re`) coordinates.
"""
mutable struct Region
    id::Int32
    cnt::Int32
    rid::Int32
    score::Int32
    qs::Int32
    qe::Int32
    rs::Int32
    re::Int32
    parent::Int32
    subsc::Int32
    as::Int32
    mlen::Int32
    blen::Int32
    n_sub::Int32
    score0::Int32
    mapq::UInt8
    split::UInt8
    rev::Bool
    inv::Bool
    sam_pri::Bool
    proper_frag::Bool
    pe_thru::Bool
    seg_split::Bool
    seg_id::UInt8
    split_inv::Bool
    is_alt::Bool
    strand_retained::Bool
    is_spliced::Bool
    hash::UInt32
    div::Float32
    extra::Union{Nothing,AlignmentExtra}
end
Region() = Region(0, 0, 0, 0, 0, 0, 0, 0, PARENT_UNSET, 0, 0, 0, 0, 0, 0, 0, 0, false, false, false, false, false,
                  false, 0, false, false, false, false, 0, -1.0f0, nothing)
Base.copy(r::Region) = Region((getfield(r, f) for f in fieldnames(Region))...)

function _fuzzy_length!(r::Region, a::Vector{MM128})
    r.mlen = r.blen = 0
    r.cnt <= 0 && return r
    first = a[r.as + 1]
    r.mlen = r.blen = _span(first)
    for i in (r.as + 1):(r.as + r.cnt - 1)
        span = _span(a[i + 1])
        tl = _ax(a[i + 1]) - _ax(a[i])
        ql = _ay(a[i + 1]) - _ay(a[i])
        r.blen += tl > ql ? tl : ql
        r.mlen += (tl > span && ql > span) ? span : (tl < ql ? tl : ql)
    end
    return r
end

# mm_reg_set_coor
function _set_coordinates!(r::Region, qlen::Int, a::Vector{MM128}, is_qstrand::Bool)
    k = r.as
    q_span = _span(a[k + 1])
    r.rev = (a[k + 1].x >> 63) == 1
    r.rid = Int32((a[k + 1].x << 1) >> 33)
    first_x = _ax(a[k + 1]); last = a[k + r.cnt]
    r.rs = first_x + Int32(1) > q_span ? first_x + Int32(1) - q_span : Int32(0)
    r.re = _ax(last) + Int32(1)
    if !r.rev || is_qstrand
        r.qs = _ay(a[k + 1]) + Int32(1) - q_span
        r.qe = _ay(last) + Int32(1)
    else
        r.qs = Int32(qlen) - (_ay(last) + Int32(1))
        r.qe = Int32(qlen) - (_ay(a[k + 1]) + Int32(1) - q_span)
    end
    return _fuzzy_length!(r, a)
end

@inline function _hash64_unmasked(key::UInt64)
    key = ~key + (key << 21)
    key = key ⊻ (key >> 24)
    key = (key + (key << 3)) + (key << 8)
    key = key ⊻ (key >> 14)
    key = (key + (key << 2)) + (key << 4)
    key = key ⊻ (key >> 28)
    key = key + (key << 31)
    return key
end

# mm_gen_regs: regions sorted by chain score, ties broken by a hash of the first
# anchor and the read.
function generate_regions(hash::UInt32, qlen::Int, u::Vector{UInt64}, a::Vector{MM128}, is_qstrand::Bool)
    n_u = length(u)
    n_u == 0 && return Region[]
    z = Vector{MM128}(undef, n_u)
    k = 0
    for i in 1:n_u
        h = UInt32(_hash64_unmasked((_hash64_unmasked(a[k + 1].x) + _hash64_unmasked(a[k + 1].y)) ⊻ UInt64(hash)) & 0xffffffff)
        z[i] = MM128(u[i] ⊻ UInt64(h), UInt64(k) << 32 | UInt64(UInt32(u[i] & 0xffffffff)))
        k += Int(u[i] & 0xffffffff)
    end
    radix_sort_128x!(z)
    reverse!(z)
    regions = Vector{Region}(undef, n_u)
    for i in 1:n_u
        r = Region()
        r.id = i - 1
        r.parent = PARENT_UNSET
        r.score = r.score0 = Int32(z[i].x >> 32)
        r.hash = UInt32(z[i].x & 0xffffffff)
        r.cnt = reinterpret(Int32, UInt32(z[i].y & 0xffffffff))
        r.as = Int32(z[i].y >> 32)
        r.div = -1.0f0
        _set_coordinates!(r, qlen, a, is_qstrand)
        regions[i] = r
    end
    return regions
end

# khash string and integer hashes used for the per-read seed.
function _x31_hash(name::AbstractString)
    bytes = codeunits(name)
    isempty(bytes) && return UInt32(0)
    h = UInt32(reinterpret(Int8, bytes[1]) % UInt32)
    h == 0 && return h
    for c in @view bytes[2:end]
        h = (h << 5) - h + (reinterpret(Int8, c) % UInt32)
    end
    return h
end

function _wang_hash(key::UInt32)
    key += ~(key << 15)
    key ⊻= key >> 10
    key += key << 3
    key ⊻= key >> 6
    key += ~(key << 11)
    key ⊻= key >> 16
    return key
end
