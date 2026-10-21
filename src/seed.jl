# Seed collection (seed.c and the seed part of map.c).

# Anchor y-field flags (mmpriv.h).
const SEED_SELF = UInt64(1) << 43
const SEED_TANDEM = UInt64(1) << 42
const SEED_SEG_SHIFT = 48
const SEED_SEG_MASK = UInt64(0xff) << 48

# mm_seed_mz_flt: drop query minimizers occurring too often in the query.
function _seed_mz_filter!(minimizers::Vector{MM128}, q_occ_max::Int, q_occ_frac::Float32)
    n = length(minimizers)
    (n <= q_occ_max || q_occ_frac <= 0 || q_occ_max <= 0) && return minimizers
    a = [MM128(minimizers[i].x, UInt64(i - 1)) for i in 1:n]
    radix_sort_128x!(a)
    keep = trues(n)
    st = 1
    for i in 2:(n + 1)
        if i == n + 1 || a[i].x != a[st].x
            count = i - st
            if count > q_occ_max && count > n * q_occ_frac
                for j in st:(i - 1)
                    keep[a[j].y + 1] = false
                end
            end
            st = i
        end
    end
    # Upstream zeroes x and then drops entries with x == 0.
    filtered = MM128[]
    for i in 1:n
        (keep[i] && minimizers[i].x != 0) && push!(filtered, minimizers[i])
    end
    resize!(minimizers, length(filtered)); copyto!(minimizers, filtered)
    return minimizers
end

mutable struct _Seed
    n::Int
    q_pos::UInt32
    q_span::UInt32
    flt::Bool
    seg_id::UInt32
    is_tandem::Bool
    cr::Vector{UInt64}
end

const _MAX_MAX_HIGH_OCC = 128

# mm_seed_select: within runs of high-occurrence seeds keep the rarest few.
function _seed_select!(seeds::Vector{_Seed}, len::Int, max_occ::Int, max_max_occ::Int, dist::Int)
    n = length(seeds)
    n <= 1 && return seeds
    any(s -> s.n > max_occ, seeds) || return seeds
    heap = Vector{UInt64}(undef, _MAX_MAX_HIGH_OCC)
    last0 = -1
    for i in 0:n
        if i == n || seeds[i + 1].n <= max_occ
            if i - last0 > 1
                ps = last0 < 0 ? 0 : Int(seeds[last0 + 1].q_pos >> 1)
                pe = i == n ? len : Int(seeds[i + 1].q_pos >> 1)
                st = last0 + 1; en = i
                max_high_occ = unsafe_trunc(Int32, Float64(pe - ps) / dist + 0.499)
                if max_high_occ > 0
                    max_high_occ = min(max_high_occ, _MAX_MAX_HIGH_OCC)
                    j = st; k = 0
                    while j < en && k < max_high_occ
                        heap[k + 1] = UInt64(seeds[j + 1].n) << 32 | UInt64(j)
                        j += 1; k += 1
                    end
                    _heapmake!(heap, k, <)
                    while j < en
                        if seeds[j + 1].n < Int32(heap[1] >> 32)
                            heap[1] = UInt64(seeds[j + 1].n) << 32 | UInt64(j)
                            _heapdown!(heap, 0, k, <)
                        end
                        j += 1
                    end
                    for j in 1:k
                        seeds[UInt32(heap[j] & 0xffffffff) + 1].flt = true
                    end
                end
                for j in st:(en - 1)
                    seeds[j + 1].flt = !seeds[j + 1].flt
                end
                for j in st:(en - 1)
                    seeds[j + 1].n > max_max_occ && (seeds[j + 1].flt = true)
                end
            end
            last0 = i
        end
    end
    return seeds
end

# mm_collect_matches: seeds of the query minimizers present in the index.
function _collect_matches(index::MinimapIndex, minimizers::Vector{MM128}, qlen::Int, max_occ::Int,
                          max_max_occ::Int, dist::Int)
    seeds = _Seed[]
    n = length(minimizers)
    for i in 1:n
        p = minimizers[i]
        positions = index_get(index, p.x >> 8)
        isempty(positions) && continue
        tandem = (i > 1 && p.x >> 8 == minimizers[i - 1].x >> 8) || (i < n && p.x >> 8 == minimizers[i + 1].x >> 8)
        push!(seeds, _Seed(length(positions), UInt32(p.y & 0xffffffff), UInt32(p.x & 0xff), false,
                           UInt32(p.y >> 32), tandem, positions))
    end
    if dist > 0 && max_max_occ > max_occ
        _seed_select!(seeds, qlen, max_occ, max_max_occ, dist)
    else
        for seed in seeds
            seed.n > max_occ && (seed.flt = true)
        end
    end
    kept = _Seed[]; mini_pos = UInt64[]
    rep_len = 0; rep_st = 0; rep_en = 0; n_a = 0
    for seed in seeds
        if seed.flt
            en = Int(seed.q_pos >> 1) + 1; st = en - Int(seed.q_span)
            if st > rep_en
                rep_len += rep_en - rep_st
                rep_st = st; rep_en = en
            else
                rep_en = en
            end
        else
            n_a += seed.n
            push!(mini_pos, UInt64(seed.q_span) << 32 | UInt64(seed.q_pos >> 1))
            push!(kept, seed)
        end
    end
    rep_len += rep_en - rep_st
    return kept, n_a, rep_len, mini_pos
end

function _skip_seed(flag::Int64, r::UInt64, seed::_Seed, qname, qlen::Int, index::MinimapIndex)
    is_self = false
    if qname !== nothing && (flag & (F_NO_DIAG | F_NO_DUAL)) != 0
        rid = Int(r >> 32)
        cmp_result = cmp(qname, index.names[rid + 1])
        if (flag & F_NO_DIAG) != 0 && cmp_result == 0 && index.lengths[rid + 1] == qlen
            UInt32(r & 0xffffffff) >> 1 == seed.q_pos >> 1 && return true, is_self
            (r & 1) == (seed.q_pos & 1) && (is_self = true)
        end
        (flag & F_NO_DUAL) != 0 && cmp_result > 0 && return true, is_self
    end
    if (flag & (F_FOR_ONLY | F_REV_ONLY)) != 0
        if (r & 1) == (seed.q_pos & 1)
            (flag & F_REV_ONLY) != 0 && return true, is_self
        else
            (flag & F_FOR_ONLY) != 0 && return true, is_self
        end
    end
    return false, is_self
end

# collect_seed_hits: anchors sorted by reference position. Anchor fields:
# x = strand << 63 | rid << 32 | reference end; y = flags | span << 32 | query end,
# both ends 0-based and inclusive, query positions on the strand of the anchor.
function collect_seed_hits(index::MinimapIndex, opt::MapOptions, max_occ::Int, minimizers::Vector{MM128},
                           qlen::Int, qname)
    seeds, n_a, rep_len, mini_pos = _collect_matches(index, minimizers, qlen, max_occ, opt.max_max_occ, opt.occ_dist)
    (opt.flag & F_HEAP_SORT) != 0 && throw(ArgumentError("heap-sorted seed collection (sr preset) is not implemented"))
    anchors = MM128[]
    sizehint!(anchors, n_a)
    for seed in seeds
        for r in seed.cr
            skip, is_self = _skip_seed(opt.flag, r, seed, qname, qlen, index)
            skip && continue
            rpos = UInt64(UInt32(r & 0xffffffff) >> 1)
            if (r & 1) == (seed.q_pos & 1)
                x = (r & 0xffffffff00000000) | rpos
                y = UInt64(seed.q_span) << 32 | UInt64(seed.q_pos >> 1)
            elseif (opt.flag & F_QSTRAND) == 0
                x = UInt64(1) << 63 | (r & 0xffffffff00000000) | rpos
                y = UInt64(seed.q_span) << 32 | UInt64(qlen - (Int(seed.q_pos >> 1) + 1 - Int(seed.q_span)) - 1)
            else
                len = index.lengths[Int(r >> 32) + 1]
                x = UInt64(1) << 63 | (r & 0xffffffff00000000) | UInt64(len - (Int(rpos) + 1 - Int(seed.q_span)) - 1)
                y = UInt64(seed.q_span) << 32 | UInt64(seed.q_pos >> 1)
            end
            y |= UInt64(seed.seg_id) << SEED_SEG_SHIFT
            seed.is_tandem && (y |= SEED_TANDEM)
            is_self && (y |= SEED_SELF)
            push!(anchors, MM128(x, y))
        end
    end
    radix_sort_128x!(anchors)
    return anchors, rep_len, mini_pos
end
