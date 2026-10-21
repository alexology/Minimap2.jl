# Primary/secondary selection, divergence estimate and mapping quality (hit.c,
# esterr.c). Region ids and parents are 0-based, as in minimap2.

_alt_score(score::Integer, alt_diff_frac::Float32) =
    score < 0 ? Int32(score) : max(unsafe_trunc(Int32, score * (1.0 - Float64(alt_diff_frac)) + 0.499), Int32(1))

# mm_set_parent: mark secondary regions and their best alternative scores.
function set_parent!(regions::Vector{Region}, mask_level::Float32, mask_len::Int, sub_diff::Int,
                     hard_mask_level::Bool, alt_diff_frac::Float32)
    n = length(regions)
    n == 0 && return regions
    for (i, r) in enumerate(regions)
        r.id = i - 1
    end
    w = Int[0]; regions[1].parent = 0
    for i in 1:(n - 1)
        ri = regions[i + 1]
        si = ri.qs; ei = ri.qe
        uncov_len = 0
        is_new_primary = false
        if !hard_mask_level
            cov = UInt64[]
            for j in w
                rp = regions[j + 1]
                sj = rp.qs; ej = rp.qe
                (ej <= si || sj >= ei) && continue
                sj < si && (sj = si); ej > ei && (ej = ei)
                push!(cov, UInt64(reinterpret(UInt32, sj)) << 32 | UInt64(reinterpret(UInt32, ej)))
            end
            if isempty(cov)
                is_new_primary = true
            else
                _radix_sort!(cov)
                x = si
                for c in cov
                    start = reinterpret(Int32, UInt32(c >> 32)); stop = reinterpret(Int32, UInt32(c & 0xffffffff))
                    start > x && (uncov_len += start - x)
                    x = stop > x ? stop : x
                end
                ei > x && (uncov_len += ei - x)
            end
        end
        found_parent = false
        if !is_new_primary
            for j in w
                rp = regions[j + 1]
                sj = rp.qs; ej = rp.qe
                (ej <= si || sj >= ei) && continue
                lmin = min(ej - sj, ei - si); lmax = max(ej - sj, ei - si)
                ol = si < sj ? (ei < sj ? 0 : ei < ej ? ei - sj : ej - sj) : (ej < si ? 0 : ej < ei ? ej - si : ei - si)
                if Float32(ol) / Float32(lmin) - Float32(uncov_len) / Float32(lmax) > mask_level && uncov_len <= mask_len
                    cnt_sub = false
                    sci = ri.score
                    ri.parent = rp.parent
                    (!rp.is_alt && ri.is_alt) && (sci = _alt_score(sci, alt_diff_frac))
                    rp.subsc = max(rp.subsc, sci)
                    ri.cnt >= rp.cnt && (cnt_sub = true)
                    if rp.extra !== nothing && ri.extra !== nothing &&
                       (rp.rid != ri.rid || rp.rs != ri.rs || rp.re != ri.re || ol != lmin)
                        sci = ri.extra.dp_max
                        (!rp.is_alt && ri.is_alt) && (sci = _alt_score(sci, alt_diff_frac))
                        rp.extra.dp_max2 = max(rp.extra.dp_max2, sci)
                        rp.extra.dp_max - ri.extra.dp_max <= sub_diff && (cnt_sub = true)
                    end
                    cnt_sub && (rp.n_sub += 1)
                    found_parent = true
                    break
                end
            end
        end
        if !found_parent
            push!(w, i); ri.parent = i; ri.n_sub = 0
        end
    end
    return regions
end

function set_sam_primary!(regions::Vector{Region})
    n_pri = 0
    for r in regions
        if r.id == r.parent
            n_pri += 1
            r.sam_pri = n_pri == 1
        else
            r.sam_pri = false
        end
    end
    return n_pri
end

# mm_sync_regs: renumber ids and remap parents after regions were removed.
function _sync_regions!(regions::Vector{Region})
    isempty(regions) && return regions
    max_id = maximum(r.id for r in regions)
    tmp = fill(-1, max_id + 1)
    for (i, r) in enumerate(regions)
        r.id >= 0 && (tmp[r.id + 1] = i - 1)
    end
    for (i, r) in enumerate(regions)
        r.id = i - 1
        if r.parent == PARENT_TMP_PRI
            r.parent = i - 1
        elseif r.parent >= 0 && tmp[r.parent + 1] >= 0
            r.parent = tmp[r.parent + 1]
        else
            r.parent = PARENT_UNSET
        end
    end
    set_sam_primary!(regions)
    return regions
end

# mm_select_sub: keep primaries and the best secondaries. Like upstream, the
# array is compacted in place while parents are read by their original index,
# which can already hold a later region.
function select_sub!(regions::Vector{Region}, pri_ratio::Float32, min_diff::Int, best_n::Int,
                     check_strand::Bool, min_strand_sc::Int)
    (pri_ratio > 0 && !isempty(regions)) || return regions
    n = length(regions)
    k = 0; n_2nd = 0
    for i in 0:(n - 1)
        ri = regions[i + 1]; p = ri.parent
        rp = regions[p + 1]
        if p == i || ri.inv
            regions[k + 1] = ri; k += 1
        elseif (ri.score >= rp.score * pri_ratio || ri.score + min_diff >= rp.score) && n_2nd < best_n
            if !(ri.qs == rp.qs && ri.qe == rp.qe && ri.rid == rp.rid && ri.rs == rp.rs && ri.re == rp.re)
                regions[k + 1] = ri; k += 1; n_2nd += 1
            end
        elseif check_strand && n_2nd < best_n && ri.score > min_strand_sc && ri.rev != rp.rev
            ri.strand_retained = true
            regions[k + 1] = ri; k += 1; n_2nd += 1
        end
    end
    if k != n
        resize!(regions, k)
        _sync_regions!(regions)
    end
    return regions
end

function filter_strand_retained!(regions::Vector{Region})
    k = 0
    for i in 0:(length(regions) - 1)
        r = regions[i + 1]
        if !r.strand_retained || r.div < regions[r.parent + 1].div * 5.0f0 || r.div < 0.01f0
            regions[k + 1] = r; k += 1
        end
    end
    resize!(regions, k)
    return regions
end

# mm_est_err: sequence divergence estimated from the fraction of query
# minimizers found in each chain.
function estimate_divergence!(regions::Vector{Region}, index::MinimapIndex, qlen::Int, a::Vector{MM128},
                              mini_pos::Vector{UInt64})
    n = length(mini_pos)
    n == 0 && return regions
    sum_k = UInt64(0)
    for m in mini_pos
        sum_k += (m >> 32) & 0xff
    end
    avg_k = Float32(sum_k) / Float32(n)
    forward_position(anchor) = (x = _ay(anchor); q_span = _span(anchor);
                                anchor.x >> 63 == 1 ? Int32(qlen) - Int32(1) - (x + Int32(1) - q_span) : x)
    for r in regions
        r.div = -1.0f0
        r.cnt == 0 && continue
        x = forward_position(r.rev ? a[r.as + r.cnt] : a[r.as + 1])
        lo = 0; hi = n - 1; st = -1
        while lo <= hi
            m = (lo + hi) >>> 1
            y = reinterpret(Int32, UInt32(mini_pos[m + 1] & 0xffffffff))
            if y < x
                lo = m + 1
            elseif y > x
                hi = m - 1
            else
                st = m
                break
            end
        end
        st < 0 && continue
        en = st
        l_ref = index.lengths[r.rid + 1]
        k = 1; j = st + 1; n_match = 1
        while j < n && k < r.cnt
            x = forward_position(r.rev ? a[r.as + r.cnt - k] : a[r.as + k + 1])
            if x == reinterpret(Int32, UInt32(mini_pos[j + 1] & 0xffffffff))
                k += 1; en = j; n_match += 1
            end
            j += 1
        end
        n_tot = en - st + 1
        (r.qs > avg_k && r.rs > avg_k) && (n_tot += 1)
        (qlen - r.qs > avg_k && l_ref - r.re > avg_k) && (n_tot += 1)
        r.div = n_match >= n_tot ? 0.0f0 : Float32(1.0 - (Float64(n_match) / n_tot)^(1.0 / Float64(avg_k)))
    end
    return regions
end

# mm_set_mapq2
function set_mapq!(regions::Vector{Region}, min_chain_sc::Int, match_sc::Int, rep_len::Int, is_sr::Bool, is_splice::Bool)
    isempty(regions) && return regions
    q_coef = 40.0f0
    sum_sc = Int64(0); n_2nd_splice = 0
    for r in regions
        if r.parent == r.id
            sum_sc += r.score
        elseif r.is_spliced
            n_2nd_splice += 1
        end
    end
    uniq_ratio = Float32(sum_sc) / Float32(sum_sc + rep_len)
    for r in regions
        if r.inv
            r.mapq = 0
        elseif r.parent == r.id
            pen_s1 = (r.score > 100 ? 1.0f0 : 0.01f0 * Float32(r.score)) * uniq_ratio
            pen_cm = r.cnt > 10 ? 1.0f0 : 0.1f0 * Float32(r.cnt)
            pen_cm = pen_s1 < pen_cm ? pen_s1 : pen_cm
            subsc = max(r.subsc, Int32(min_chain_sc))
            p = r.extra
            if p !== nothing && p.dp_max2 > 0 && p.dp_max > 0
                identity = Float32(r.mlen) / Float32(r.blen)
                x = if is_sr && is_splice
                    Float32(p.dp_max2) / Float32(p.dp_max)
                else
                    Float32(p.dp_max2) * Float32(subsc) / Float32(p.dp_max) / Float32(r.score0)
                end
                mapq = unsafe_trunc(Int32, identity * pen_cm * q_coef * (1.0f0 - x * x) * log(Float32(p.dp_max) / Float32(match_sc)))
                if !is_sr
                    mapq_alt = unsafe_trunc(Int32, 6.02f0 * identity * identity * Float32(p.dp_max - p.dp_max2) / Float32(match_sc) + 0.499f0)
                    mapq = min(mapq, mapq_alt)
                end
                (is_splice && is_sr && r.is_spliced && n_2nd_splice == 0) && (mapq += 10)
            else
                x = Float32(subsc) / Float32(r.score0)
                if p !== nothing
                    identity = Float32(r.mlen) / Float32(r.blen)
                    mapq = unsafe_trunc(Int32, identity * pen_cm * q_coef * (1.0f0 - x) * log(Float32(p.dp_max) / Float32(match_sc)))
                else
                    mapq = unsafe_trunc(Int32, pen_cm * q_coef * (1.0f0 - x) * log(Float32(r.score)))
                end
            end
            mapq -= unsafe_trunc(Int32, 4.343f0 * log(Float32(r.n_sub + 1)) + 0.499f0)
            mapq = max(mapq, Int32(0))
            r.mapq = UInt8(min(mapq, Int32(60)))
            (p !== nothing && p.dp_max > p.dp_max2 && r.mapq == 0) && (r.mapq = 1)
        else
            r.mapq = 0
        end
    end
    _set_inversion_mapq!(regions)
    return regions
end

function _set_inversion_mapq!(regions::Vector{Region})
    n = length(regions)
    n < 3 && return regions
    any(r -> r.inv, regions) || return regions
    aux = MM128[]
    for (i, r) in enumerate(regions)
        (r.parent == i - 1 || r.parent < 0) && push!(aux, MM128(UInt64(r.rid) << 32 | UInt64(reinterpret(UInt32, r.rs)), UInt64(i - 1)))
    end
    radix_sort_128x!(aux)
    for i in 2:(length(aux) - 1)
        inv = regions[aux[i].y + 1]
        if inv.inv
            inv.mapq = min(regions[aux[i - 1].y + 1].mapq, regions[aux[i + 1].y + 1].mapq)
        end
    end
    return regions
end

# mm_hit_sort: drop soft-deleted regions and sort by DP (or chaining) score.
function hit_sort!(regions::Vector{Region}, alt_diff_frac::Float32)
    n = length(regions)
    n <= 1 && return regions
    aux = MM128[]
    for (i, r) in enumerate(regions)
        if r.inv || r.cnt > 0
            score = r.extra !== nothing ? r.extra.dp_max : r.score
            r.is_alt && (score = _alt_score(score, alt_diff_frac))
            push!(aux, MM128(UInt64(reinterpret(UInt32, Int32(score))) << 32 | UInt64(r.hash), UInt64(i - 1)))
        end
    end
    radix_sort_128x!(aux)
    sorted = [regions[aux[i].y + 1] for i in length(aux):-1:1]
    resize!(regions, length(sorted)); copyto!(regions, sorted)
    return regions
end

# mm_filter_regs (only meaningful after base-level alignment).
function filter_regions!(regions::Vector{Region}, opt::MapOptions, qlen::Int)
    kept = Region[]
    for r in regions
        flt = !r.inv && !r.seg_split && r.cnt < opt.min_cnt
        if r.extra !== nothing
            if r.mlen < opt.min_chain_score
                flt = true
            elseif r.extra.dp_max < opt.min_dp_max
                flt = true
            elseif r.qs > qlen * opt.max_clip_ratio && qlen - r.qe > qlen * opt.max_clip_ratio
                flt = true
            end
        end
        flt || push!(kept, r)
    end
    resize!(regions, length(kept)); copyto!(regions, kept)
    return regions
end
