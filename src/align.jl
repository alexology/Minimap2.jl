# Base-level alignment of chains (align.c, `mm_align_skeleton`), non-spliced.

const F_EQX = Int64(0x4000000)
const SEED_IGNORE = UInt64(1) << 41
const SEED_LONG_JOIN = UInt64(1) << 40

# ksw_gen_ts_mat for m = 5 (A, C, G, T, N).
function _score_matrix(a::Int, b::Int, transition::Int, sc_ambi::Int)
    m = 5
    a = abs(a); b = b > 0 ? -b : b; sc_ambi = sc_ambi > 0 ? -sc_ambi : sc_ambi
    mat = Vector{Int8}(undef, 25)
    for i in 0:(m - 2)
        for j in 0:(m - 2)
            mat[i * m + j + 1] = i == j ? a : b
        end
        mat[i * m + m] = sc_ambi
    end
    for j in 0:(m - 1)
        mat[(m - 1) * m + j + 1] = sc_ambi
    end
    if transition != 0 && transition != b
        transition = transition > 0 ? -transition : transition
        mat[0 * m + 2 + 1] = transition; mat[1 * m + 3 + 1] = transition
        mat[2 * m + 0 + 1] = transition; mat[3 * m + 1 + 1] = transition
    end
    return mat
end

function _align_pair!(ez::KswResult, opt::MapOptions, qseq::AbstractVector{UInt8}, tseq::AbstractVector{UInt8},
                      mat::Vector{Int8}, w::Int, end_bonus::Int, zdrop::Int, ksw_flag::Int)
    qlen = length(qseq); tlen = length(tseq)
    (opt.transition != 0 && opt.b != opt.transition) && (ksw_flag |= KSW_EZ_GENERIC_SC)
    if opt.max_sw_mat > 0 && Int64(tlen) * qlen > opt.max_sw_mat
        _reset!(ez); ez.zdropped = true
    elseif (opt.flag & F_SPLICE) != 0
        throw(ArgumentError("spliced alignment is not implemented"))
    elseif opt.q == opt.q2 && opt.e == opt.e2
        throw(ArgumentError("single affine gap alignment (ksw_extz2) is not implemented"))
    else
        ksw_extd2!(ez, qseq, tseq, 5, mat, opt.q, opt.e, opt.q2, opt.e2, w, zdrop, end_bonus, ksw_flag)
    end
    return ez
end

# mm_test_zdrop: 0 = fine, 1 = Z-drop, 2 = possible inversion.
function _test_zdrop(opt::MapOptions, qseq::AbstractVector{UInt8}, tseq::AbstractVector{UInt8},
                     cigar::Vector{UInt32}, mat::Vector{Int8})
    score = 0; max_score = Int(typemin(Int32)); max_i = -1; max_j = -1; i = 0; j = 0; max_zdrop = 0
    pos11 = pos12 = pos21 = pos22 = -1    # pos[0][0], pos[0][1], pos[1][0], pos[1][1]
    for c in cigar
        op = Int(c & 0xf); len = Int(c >> 4)
        if op == CIGAR_MATCH
            for l in 0:(len - 1)
                score += mat[Int(tseq[i + l + 1]) * 5 + Int(qseq[j + l + 1]) + 1]
                ii = i + l; jj = j + l
                if score < max_score
                    li = ii - max_i; lj = jj - max_j
                    z = max_score - score - abs(li - lj) * opt.e
                    if z > max_zdrop
                        max_zdrop = z; pos11 = max_i; pos12 = ii; pos21 = max_j; pos22 = jj
                    end
                else
                    max_score = score; max_i = ii; max_j = jj
                end
            end
            i += len; j += len
        elseif op == CIGAR_INS || op == CIGAR_DEL || op == CIGAR_N_SKIP
            score -= opt.q + opt.e * len
            op == CIGAR_INS ? (j += len) : (i += len)
            if score < max_score
                li = i - max_i; lj = j - max_j
                z = max_score - score - abs(li - lj) * opt.e
                if z > max_zdrop
                    max_zdrop = z; pos11 = max_i; pos12 = i; pos21 = max_j; pos22 = j
                end
            else
                max_score = score; max_i = i; max_j = j
            end
        end
    end
    q_len = pos22 - pos21; t_len = pos12 - pos11
    if (opt.flag & (F_SPLICE | F_SR | F_FOR_ONLY | F_REV_ONLY)) == 0 && max_zdrop > opt.zdrop_inv &&
       q_len < opt.max_gap && t_len < opt.max_gap
        qseq2 = UInt8[(c = qseq[pos22 - k]; c >= 4 ? 0x04 : 0x03 - c) for k in 0:(q_len - 1)]
        inv_score, _, _ = ksw_ll_i16(qseq2, view(tseq, (pos11 + 1):(pos11 + t_len)), 5, mat, opt.q, opt.e)
        (inv_score >= opt.min_chain_score * opt.a && inv_score >= opt.min_dp_max) && return 2
    end
    return max_zdrop > opt.zdrop ? 1 : 0
end

_new_extra() = AlignmentExtra(0, 0, 0, 0, 0, 0, UInt32[])

function _append_cigar!(r::Region, cigar::Vector{UInt32})
    isempty(cigar) && return r
    r.extra === nothing && (r.extra = _new_extra())
    p = r.extra
    if !isempty(p.cigar) && (p.cigar[end] & 0xf) == (cigar[1] & 0xf)
        p.cigar[end] += (cigar[1] >> 4) << 4
        append!(p.cigar, @view cigar[2:end])
    else
        append!(p.cigar, cigar)
    end
    return r
end

# mm_fix_cigar: left-align indels, merge adjacent I/D runs, drop a leading gap.
function _fix_cigar!(r::Region, qseq::AbstractVector{UInt8}, tseq::AbstractVector{UInt8})
    p = r.extra; cigar = p.cigar
    qshift = tshift = 0
    length(cigar) <= 1 && return qshift, tshift
    toff = 0; qoff = 0; to_shrink = false
    n = length(cigar)
    for k in 1:n
        op = Int(cigar[k] & 0xf); len = Int(cigar[k] >> 4)
        len == 0 && (to_shrink = true)
        if op == CIGAR_MATCH
            toff += len; qoff += len
        elseif op == CIGAR_INS || op == CIGAR_DEL
            if k > 1 && k < n && (cigar[k - 1] & 0xf) == 0 && (cigar[k + 1] & 0xf) == 0
                prev_len = Int(cigar[k - 1] >> 4)
                l = 0
                if op == CIGAR_INS
                    while l < prev_len && qseq[qoff - 1 - l + 1] == qseq[qoff + len - 1 - l + 1]
                        l += 1
                    end
                else
                    while l < prev_len && tseq[toff - 1 - l + 1] == tseq[toff + len - 1 - l + 1]
                        l += 1
                    end
                end
                if l > 0
                    cigar[k - 1] -= UInt32(l) << 4; cigar[k + 1] += UInt32(l) << 4; qoff -= l; toff -= l
                end
                l == prev_len && (to_shrink = true)
            end
            op == CIGAR_INS ? (qoff += len) : (toff += len)
        elseif op == CIGAR_N_SKIP
            toff += len
        end
    end
    @assert qoff == r.qe - r.qs && toff == r.re - r.rs
    k = 1
    while k <= n - 2
        if (cigar[k] & 0xf) > 0 && (cigar[k] & 0xf) + (cigar[k + 1] & 0xf) == 3
            s = [0, 0, 0]
            l = k
            while l <= n
                op = Int(cigar[l] & 0xf)
                if op == CIGAR_INS || op == CIGAR_DEL || (cigar[l] >> 4) == 0
                    s[op + 1] += Int(cigar[l] >> 4)
                else
                    break
                end
                l += 1
            end
            if s[2] > 0 && s[3] > 0 && l - k > 2
                cigar[k] = UInt32(s[2]) << 4 | UInt32(CIGAR_INS)
                cigar[k + 1] = UInt32(s[3]) << 4 | UInt32(CIGAR_DEL)
                kk = k + 2
                while kk < l
                    cigar[kk] &= 0xf; kk += 1
                end
                to_shrink = true
            end
            k = l
        end
        k += 1
    end
    if to_shrink
        filter!(c -> (c >> 4) != 0, cigar)
        merged = UInt32[]
        for k in eachindex(cigar)
            if k == length(cigar) || (cigar[k] & 0xf) != (cigar[k + 1] & 0xf)
                push!(merged, cigar[k])
            else
                cigar[k + 1] += (cigar[k] >> 4) << 4
            end
        end
        empty!(cigar); append!(cigar, merged)
    end
    op0 = Int(cigar[1] & 0xf)
    if op0 == CIGAR_INS || op0 == CIGAR_DEL
        l = Int(cigar[1] >> 4)
        if op0 == CIGAR_INS
            r.rev ? (r.qe -= l) : (r.qs += l)
            qshift = l
        else
            r.rs += l; tshift = l
        end
        deleteat!(cigar, 1)
    end
    return qshift, tshift
end

function _update_cigar_eqx!(r::Region, qseq::AbstractVector{UInt8}, tseq::AbstractVector{UInt8})
    p = r.extra
    p === nothing && return r
    out = UInt32[]; toff = 0; qoff = 0
    for c in p.cigar
        op = Int(c & 0xf); len = Int(c >> 4)
        if op == CIGAR_MATCH
            while len > 0
                l = 0
                while l < len && qseq[qoff + l + 1] == tseq[toff + l + 1]
                    l += 1
                end
                l > 0 && push!(out, UInt32(l) << 4 | 0x7)
                len -= l; toff += l; qoff += l
                l = 0
                while l < len && qseq[qoff + l + 1] != tseq[toff + l + 1]
                    l += 1
                end
                l > 0 && push!(out, UInt32(l) << 4 | 0x8)
                len -= l; toff += l; qoff += l
            end
            continue
        elseif op == CIGAR_INS
            qoff += len
        elseif op == CIGAR_DEL || op == CIGAR_N_SKIP
            toff += len
        end
        push!(out, c)
    end
    p.cigar = out
    return r
end

# mm_update_extra: fix the CIGAR and recompute matches, block length and DP max.
function _update_extra!(r::Region, qseq::AbstractVector{UInt8}, tseq::AbstractVector{UInt8}, mat::Vector{Int8},
                        q::Int, e::Int, is_eqx::Bool, log_gap::Bool)
    p = r.extra
    p === nothing && return r
    qshift, tshift = _fix_cigar!(r, qseq, tseq)
    qseq = view(qseq, (qshift + 1):length(qseq)); tseq = view(tseq, (tshift + 1):length(tseq))
    r.blen = r.mlen = 0; r.is_spliced = false
    s = 0.0; best = 0.0; toff = 0; qoff = 0
    for c in p.cigar
        op = Int(c & 0xf); len = Int(c >> 4)
        if op == CIGAR_MATCH
            n_ambi = 0; n_diff = 0
            for l in 0:(len - 1)
                cq = Int(qseq[qoff + l + 1]); ct = Int(tseq[toff + l + 1])
                if ct > 3 || cq > 3
                    n_ambi += 1
                elseif ct != cq
                    n_diff += 1
                end
                s += mat[ct * 5 + cq + 1]
                if s < 0
                    s = 0.0
                else
                    best = best > s ? best : s
                end
            end
            r.blen += len - n_ambi; r.mlen += len - (n_ambi + n_diff); p.n_ambi += n_ambi
            toff += len; qoff += len
        elseif op == CIGAR_INS || op == CIGAR_DEL
            n_ambi = 0
            for l in 0:(len - 1)
                code = op == CIGAR_INS ? qseq[qoff + l + 1] : tseq[toff + l + 1]
                code > 3 && (n_ambi += 1)
            end
            r.blen += len - n_ambi; p.n_ambi += n_ambi
            s -= log_gap ? q + Float64(e) * Float64(_mg_log2(Float32(1.0 + len))) : q + e
            s < 0 && (s = 0.0)
            op == CIGAR_INS ? (qoff += len) : (toff += len)
        elseif op == CIGAR_N_SKIP
            r.is_spliced = true; toff += len
        end
    end
    p.dp_max = p.dp_max0 = unsafe_trunc(Int32, best + 0.499)
    @assert qoff == r.qe - r.qs && toff == r.re - r.rs
    is_eqx && _update_cigar_eqx!(r, qseq, tseq)
    return r
end

# mm_get_hplen_back
function _hplen_back(index::MinimapIndex, rid::Int, x::Int)
    off0 = index.offsets[rid + 1]; off = off0 + x
    c = index.codes[off + 1]
    i = off - 1
    while i >= off0
        index.codes[i + 1] != c && break
        i -= 1
    end
    return off - i
end

# mm_adjust_minier: start position of the minimizer of an anchor.
function _adjust_minimizer(index::MinimapIndex, qseqs, a::MM128)
    if index.hpc
        qseq = qseqs[Int(a.x >> 63) + 1]
        q = Int(_ay(a))
        c = qseq[q + 1]
        i = q - 1
        while i > 0
            qseq[i + 1] != c && break
            i -= 1
        end
        q = i + 1
        c2 = _hplen_back(index, Int((a.x << 1) >> 33), Int(_ax(a)))
        return Int(_ax(a)) + 1 - c2, q
    end
    return Int(_ax(a)) - (index.k >> 1), Int(_ay(a)) - (index.k >> 1)
end

_gap_of(a, i) = (Int(_ay(a[i])) - Int(_ay(a[i - 1]))) - (Int(_ax(a[i])) - Int(_ax(a[i - 1])))
_gap_of_u(a, i) = (Int(_ay(a[i])) - Int(_ay(a[i - 1]))) - Int(reinterpret(Int32, UInt32((a[i].x - a[i - 1].x) & 0xffffffff)))

function _collect_long_gaps(as1::Int, cnt1::Int, a::Vector{MM128}, min_gap::Int)
    K = Int[]
    for i in 1:(cnt1 - 1)
        gap = _gap_of(a, as1 + i + 1)
        (gap < -min_gap || gap > min_gap) && push!(K, i)
    end
    return length(K) <= 1 ? Int[] : K
end

function _set_flag!(a::Vector{MM128}, i::Int, flag::UInt64)
    a[i] = MM128(a[i].x, a[i].y | flag)
end

function _filter_bad_seeds!(as1::Int, cnt1::Int, a::Vector{MM128}, min_gap::Int, diff_thres::Int,
                            max_ext_len::Int, max_ext_cnt::Int)
    K = _collect_long_gaps(as1, cnt1, a, min_gap)
    isempty(K) && return a
    n = length(K)
    best = 0; max_st = max_en = -1
    k = 0
    while true
        if k == n || k >= max_en
            if max_en > 0
                for i in K[max_st + 1]:(K[max_en + 1] - 1)
                    _set_flag!(a, as1 + i + 1, SEED_IGNORE)
                end
            end
            best = 0; max_st = max_en = -1
            k == n && break
        end
        i = K[k + 1]
        n_ins = 0; n_del = 0
        gap = _gap_of_u(a, as1 + i + 1)
        gap > 0 ? (n_ins += gap) : (n_del += -gap)
        qs = Int(_ay(a[as1 + i])); rs = Int(_ax(a[as1 + i]))
        max_diff = 0; max_diff_l = -1
        l = k + 1
        while l < n && l <= k + max_ext_cnt
            j = K[l + 1]
            (Int(_ay(a[as1 + j + 1])) - qs > max_ext_len || Int(_ax(a[as1 + j + 1])) - rs > max_ext_len) && break
            gap = _gap_of_u(a, as1 + j + 1)
            gap > 0 ? (n_ins += gap) : (n_del += -gap)
            diff = n_ins + n_del - abs(n_ins - n_del)
            if max_diff < diff
                max_diff = diff; max_diff_l = l
            end
            l += 1
        end
        if max_diff > diff_thres && max_diff > best
            best = max_diff; max_st = k; max_en = max_diff_l
        end
        k += 1
    end
    return a
end

function _filter_bad_seeds_alt!(as1::Int, cnt1::Int, a::Vector{MM128}, min_gap::Int, max_ext::Int)
    K = _collect_long_gaps(as1, cnt1, a, min_gap)
    isempty(K) && return a
    n = length(K)
    k = 0
    while k < n
        i = K[k + 1]
        gap1 = abs(_gap_of(a, as1 + i + 1))
        re1 = Int(_ax(a[as1 + i + 1])); qe1 = Int(_ay(a[as1 + i + 1]))
        l = k + 1
        while l < n
            j = K[l + 1]
            (Int(_ay(a[as1 + j + 1])) - qe1 > max_ext || Int(_ax(a[as1 + j + 1])) - re1 > max_ext) && break
            gap2 = abs(_gap_of_u(a, as1 + j + 1))
            q_span_pre = Int(_span(a[as1 + j]))
            rs2 = Int(_ax(a[as1 + j])) + q_span_pre; qs2 = Int(_ay(a[as1 + j])) + q_span_pre
            m = min(rs2 - re1, qs2 - qe1)
            m > gap1 + gap2 && break
            re1 = Int(_ax(a[as1 + j + 1])); qe1 = Int(_ay(a[as1 + j + 1]))
            gap1 = gap2
            l += 1
        end
        if l > k + 1
            stop = K[l]
            for j in K[k + 1]:(stop - 1)
                _set_flag!(a, as1 + j + 1, SEED_IGNORE)
            end
            _set_flag!(a, as1 + stop + 1, SEED_LONG_JOIN)
        end
        k = l
    end
    return a
end

# mm_fix_bad_ends: drop poorly supported anchors at chain ends; returns (as, cnt).
function _fix_bad_ends(r::Region, a::Vector{MM128}, bw::Int, min_match::Int)
    as1 = Int(r.as); cnt = Int(r.cnt)
    r.cnt < 3 && return as1, cnt
    l = m = Int(_span(a[r.as + 1]))
    for i in (r.as + 1):(r.as + r.cnt - 2)
        q_span = Int(_span(a[i + 1]))
        (a[i + 1].y & SEED_LONG_JOIN) != 0 && break
        lr = Int(_ax(a[i + 1])) - Int(_ax(a[i])); lq = Int(_ay(a[i + 1])) - Int(_ay(a[i]))
        lo = min(lr, lq); hi = max(lr, lq)
        hi - lo > l >> 1 && (as1 = i)
        l += lo
        m += min(lo, q_span)
        (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r.mlen >> 1) && break
    end
    cnt = r.as + r.cnt - as1
    l = m = Int(_span(a[r.as + r.cnt]))
    i = r.as + r.cnt - 2
    while i > as1
        q_span = Int(_span(a[i + 2]))
        (a[i + 2].y & SEED_LONG_JOIN) != 0 && break
        lr = Int(_ax(a[i + 2])) - Int(_ax(a[i + 1])); lq = Int(_ay(a[i + 2])) - Int(_ay(a[i + 1]))
        lo = min(lr, lq); hi = max(lr, lq)
        hi - lo > l >> 1 && (cnt = i + 1 - as1)
        l += lo
        m += min(lo, q_span)
        (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r.mlen >> 1) && break
        i -= 1
    end
    return as1, cnt
end

# mm_split_reg
function _split_region!(r::Region, n::Int, qlen::Int, a::Vector{MM128}, is_qstrand::Bool)
    (n <= 0 || n >= r.cnt) && return nothing
    r2 = copy(r)
    r2.id = -1; r2.sam_pri = false; r2.extra = nothing; r2.split_inv = false
    r2.cnt = r.cnt - n
    r2.score = unsafe_trunc(Int32, Float64(Float32(r.score) * (Float32(r2.cnt) / Float32(r.cnt))) + 0.499)
    r2.as = r.as + n
    r.parent == r.id && (r2.parent = PARENT_TMP_PRI)
    _set_coordinates!(r2, qlen, a, is_qstrand)
    r.cnt -= r2.cnt; r.score -= r2.score
    _set_coordinates!(r, qlen, a, is_qstrand)
    r.split |= 1; r2.split |= 2
    return r2
end

# mm_align1: extend and fill one chain; returns the region split off by Z-drop, if any.
function _align1!(opt::MapOptions, index::MinimapIndex, qlen::Int, qseqs, r::Region, a::Vector{MM128}, ez::KswResult)
    is_sr = (opt.flag & F_SR) != 0
    r.cnt == 0 && return nothing
    rid = Int((a[r.as + 1].x << 1) >> 33); rev = (a[r.as + 1].x >> 63) == 1
    mat = _score_matrix(opt.a, opt.b, opt.transition, opt.sc_ambi)
    bw = unsafe_trunc(Int, opt.bw * 1.5 + 1.0)
    bw_long = max(unsafe_trunc(Int, opt.bw_long * 1.5 + 1.0), bw)
    is_sr && throw(ArgumentError("short-read alignment is not implemented"))
    (opt.flag & F_QSTRAND) != 0 && throw(ArgumentError("query-strand mode is not implemented"))
    if (opt.flag & F_NO_END_FLT) == 0
        as1, cnt1 = _fix_bad_ends(r, a, opt.bw, opt.min_chain_score * 2)
    else
        as1 = Int(r.as); cnt1 = Int(r.cnt)
    end
    _filter_bad_seeds!(as1, cnt1, a, 10, 40, opt.max_gap >> 1, 10)
    _filter_bad_seeds_alt!(as1, cnt1, a, 30, opt.max_gap >> 1)
    rs, qs = _adjust_minimizer(index, qseqs, a[as1 + 1])
    re, qe = _adjust_minimizer(index, qseqs, a[as1 + cnt1])
    @assert cnt1 > 0
    ksw_flag = 0
    reflen = index.lengths[rid + 1]
    # Regions for DP (start).
    first = a[r.as + 1]
    rs0 = Int(_ax(first)) + 1 - Int(_span(first)); qs0 = Int(_ay(first)) + 1 - Int(_span(first))
    rs0 < 0 && (rs0 = 0)
    rs1 = qs1 = 0
    i = r.as - 1; l = 0
    while i >= 0 && (a[i + 1].x >> 32) == (first.x >> 32)
        x = Int(_ax(a[i + 1])) + 1 - Int(_span(a[i + 1])); y = Int(_ay(a[i + 1])) + 1 - Int(_span(a[i + 1]))
        if x < rs0 && y < qs0
            l += 1
            if l > opt.min_cnt
                l = max(rs0 - x, qs0 - y)
                rs1 = rs0 - l; qs1 = qs0 - l
                rs1 < 0 && (rs1 = 0)
                break
            end
        end
        i -= 1
    end
    if qs > 0 && rs > 0
        l = min(qs, opt.max_gap)
        qs1 = max(qs1, qs - l)
        qs0 = min(qs0, qs1)
        l += l * opt.a > opt.q ? (l * opt.a - opt.q) ÷ opt.e : 0
        l = min(l, opt.max_gap); l = min(l, rs)
        rs1 = max(rs1, rs - l)
        rs0 = min(rs0, rs1); rs0 = min(rs0, rs)
    else
        rs0 = rs; qs0 = qs
    end
    # Regions for DP (end).
    last = a[r.as + r.cnt]
    re0 = Int(_ax(last)) + 1; qe0 = Int(_ay(last)) + 1
    re1 = reflen; qe1 = qlen
    i = r.as + r.cnt; l = 0
    while i < length(a) && (a[i + 1].x >> 32) == (first.x >> 32)
        x = Int(_ax(a[i + 1])) + 1; y = Int(_ay(a[i + 1])) + 1
        if x > re0 && y > qe0
            l += 1
            if l > opt.min_cnt
                l = max(x - re0, y - qe0)
                re1 = re0 + l; qe1 = qe0 + l
                break
            end
        end
        i += 1
    end
    if qe < qlen && re < reflen
        l = min(qlen - qe, opt.max_gap)
        qe1 = min(qe1, qe + l)
        qe0 = max(qe0, qe1)
        l += l * opt.a > opt.q ? (l * opt.a - opt.q) ÷ opt.e : 0
        l = min(l, opt.max_gap); l = min(l, reflen - re)
        re1 = min(re1, re + l)
        re0 = max(re0, re1)
    else
        re0 = re; qe0 = qe
    end
    if (first.y & SEED_SELF) != 0
        max_ext = abs(r.qs - r.rs)
        r.rs - rs0 > max_ext && (rs0 = r.rs - max_ext)
        r.qs - qs0 > max_ext && (qs0 = r.qs - max_ext)
        max_ext = abs(r.qe - r.re)
        re0 - r.re > max_ext && (re0 = r.re + max_ext)
        qe0 - r.qe > max_ext && (qe0 = r.qe + max_ext)
    end
    @assert re0 > rs0
    qseq_rev = qseqs[rev + 1]
    dropped = false
    r2 = nothing
    if qs > 0 && rs > 0
        qseq = reverse(qseq_rev[(qs0 + 1):qs])
        tseq = reverse(reference_codes(index, rid, rs0, rs))
        _align_pair!(ez, opt, qseq, tseq, mat, bw, opt.end_bonus, r.split_inv ? opt.zdrop_inv : opt.zdrop,
                     ksw_flag | KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR)
        if !isempty(ez.cigar)
            _append_cigar!(r, ez.cigar)
            r.extra.dp_score += ez.max
        end
        rs1 = rs - (ez.reach_end ? ez.mqe_t + 1 : ez.max_t + 1)
        qs1 = qs - (ez.reach_end ? qs - qs0 : ez.max_q + 1)
    else
        rs1 = rs; qs1 = qs
    end
    re1 = rs; qe1 = qs
    @assert qs1 >= 0 && rs1 >= 0
    for i in 1:(cnt1 - 1)
        ai = a[as1 + i + 1]
        ((ai.y & (SEED_IGNORE | SEED_TANDEM)) != 0 && i != cnt1 - 1) && continue
        re, qe = _adjust_minimizer(index, qseqs, ai)
        re1 = re; qe1 = qe
        if i == cnt1 - 1 || (ai.y & SEED_LONG_JOIN) != 0 || (qe - qs >= opt.min_ksw_len && re - rs >= opt.min_ksw_len)
            bw1 = bw_long
            (ai.y & SEED_LONG_JOIN) != 0 && (bw1 = max(qe - qs, re - rs))
            qseq = qseq_rev[(qs + 1):qe]
            tseq = reference_codes(index, rid, rs, re)
            _align_pair!(ez, opt, qseq, tseq, mat, bw1, -1, opt.zdrop, ksw_flag | KSW_EZ_APPROX_MAX)
            zdrop_code = _test_zdrop(opt, qseq, tseq, ez.cigar, mat)
            if zdrop_code != 0
                _align_pair!(ez, opt, qseq, tseq, mat, bw1, -1, zdrop_code == 2 ? opt.zdrop_inv : opt.zdrop, ksw_flag)
            end
            isempty(ez.cigar) || _append_cigar!(r, ez.cigar)
            if ez.zdropped
                r.extra === nothing && (r.extra = _new_extra())
                j = i - 1
                while j >= 0
                    Int(_ax(a[as1 + j + 1])) <= rs + ez.max_t && break
                    j -= 1
                end
                dropped = true
                j < 0 && (j = 0)
                r.extra.dp_score += ez.max
                re1 = rs + (ez.max_t + 1); qe1 = qs + (ez.max_q + 1)
                if cnt1 - (j + 1) >= opt.min_cnt
                    r2 = _split_region!(r, as1 + j + 1 - r.as, qlen, a, false)
                    (r2 !== nothing && zdrop_code == 2) && (r2.split_inv = true)
                end
                break
            else
                r.extra.dp_score += ez.score
            end
            rs = re; qs = qe
        end
    end
    if !dropped && qe < qe0 && re < re0
        qseq = qseq_rev[(qe + 1):qe0]
        tseq = reference_codes(index, rid, re, re0)
        _align_pair!(ez, opt, qseq, tseq, mat, bw, opt.end_bonus, opt.zdrop, ksw_flag | KSW_EZ_EXTZ_ONLY)
        if !isempty(ez.cigar)
            _append_cigar!(r, ez.cigar)
            r.extra.dp_score += ez.max
        end
        re1 = re + (ez.reach_end ? ez.mqe_t + 1 : ez.max_t + 1)
        qe1 = qe + (ez.reach_end ? qe0 - qe : ez.max_q + 1)
    end
    @assert qe1 <= qlen
    r.rs = rs1; r.re = re1
    if !rev
        r.qs = qs1; r.qe = qe1
    else
        r.qs = qlen - qe1; r.qe = qlen - qs1
    end
    @assert re1 - rs1 <= re0 - rs0
    if r.extra !== nothing
        tseq = reference_codes(index, rid, rs1, re1)
        qseq = qseqs[r.rev + 1][(qs1 + 1):end]
        _update_extra!(r, qseq, tseq, mat, opt.q, opt.e, (opt.flag & F_EQX) != 0, !is_sr)
    end
    return r2
end

# mm_align1_inv: align the gap between two split parts on the opposite strand.
function _align1_inv(opt::MapOptions, index::MinimapIndex, qlen::Int, qseqs, r1::Region, r2::Region, ez::KswResult)
    ((r1.split & 1) == 0 || (r2.split & 2) == 0) && return nothing
    (r1.id != r1.parent && r1.parent != PARENT_TMP_PRI) && return nothing
    (r2.id != r2.parent && r2.parent != PARENT_TMP_PRI) && return nothing
    (r1.rid != r2.rid || r1.rev != r2.rev) && return nothing
    ql = r1.rev ? r1.qs - r2.qe : r2.qs - r1.qe
    tl = r2.rs - r1.re
    (ql < opt.min_chain_score || ql > opt.max_gap) && return nothing
    (tl < opt.min_chain_score || tl > opt.max_gap) && return nothing
    mat = _score_matrix(opt.a, opt.b, opt.transition, opt.sc_ambi)
    tseq = reference_codes(index, Int(r1.rid), r1.re, r2.rs)
    qseq = r1.rev ? qseqs[1][(r2.qe + 1):(r2.qe + ql)] : qseqs[2][(qlen - r2.qs + 1):(qlen - r2.qs + ql)]
    score, q_off, t_off = ksw_ll_i16(reverse(qseq), reverse(tseq), 5, mat, opt.q, opt.e)
    score < opt.min_dp_max && return nothing
    q_off = ql - (q_off + 1); t_off = tl - (t_off + 1)
    _align_pair!(ez, opt, qseq[(q_off + 1):end], tseq[(t_off + 1):end], mat, unsafe_trunc(Int, opt.bw * 1.5), -1,
                 opt.zdrop, KSW_EZ_EXTZ_ONLY)
    isempty(ez.cigar) && return nothing
    r_inv = Region()
    r_inv.parent = 0; r_inv.div = 0.0f0
    _append_cigar!(r_inv, ez.cigar)
    r_inv.extra.dp_score = ez.max
    r_inv.id = -1
    r_inv.parent = PARENT_UNSET
    r_inv.inv = true
    r_inv.rev = !r1.rev
    r_inv.rid = r1.rid
    r_inv.div = -1.0f0
    if !r_inv.rev
        r_inv.qs = r2.qe + q_off
        r_inv.qe = r_inv.qs + ez.max_q + 1
    else
        r_inv.qe = r2.qs - q_off
        r_inv.qs = r_inv.qe - (ez.max_q + 1)
    end
    r_inv.rs = r1.re + t_off
    r_inv.re = r_inv.rs + ez.max_t + 1
    _update_extra!(r_inv, qseq[(q_off + 1):end], tseq[(t_off + 1):end], mat, opt.q, opt.e, (opt.flag & F_EQX) != 0,
                   (opt.flag & (F_SR | F_SR_RNA)) == 0)
    return r_inv
end

function _count_gaps(r::Region)
    r.extra === nothing && return -1, -1
    n_gap = 0; n_gapo = 0
    for c in r.extra.cigar
        op = Int(c & 0xf)
        if op == CIGAR_INS || op == CIGAR_DEL
            n_gapo += 1; n_gap += Int(c >> 4)
        end
    end
    return n_gap, n_gapo
end

"""
    event_identity(region) -> Float64

Gap-compressed identity of an aligned region (minimap2 `mm_event_identity`).
"""
function event_identity(r::Region)
    r.extra === nothing && return -1.0
    n_gap, n_gapo = _count_gaps(r)
    return Float64(r.mlen) / (r.blen + Int(r.extra.n_ambi) - n_gap + n_gapo)
end

function _recal_max_dp(r::Region, b2::Float64, match_sc::Int)
    r.extra === nothing && return -1
    n_gap = 0; gap_cost = 0.0
    for c in r.extra.cigar
        op = Int(c & 0xf); len = Int(c >> 4)
        if op == CIGAR_INS || op == CIGAR_DEL
            gap_cost += b2 + Float64(_mg_log2(Float32(1.0 + len)))
            n_gap += len
        end
    end
    n_mis = r.blen + Int(r.extra.n_ambi) - r.mlen - n_gap
    return unsafe_trunc(Int32, match_sc * (r.mlen - b2 * n_mis - gap_cost) + 0.499)
end

# mm_update_dp_max: rescore by gap-compressed identity when two hits are close.
function _update_dp_max!(qlen::Int, regions::Vector{Region}, frac::Float32, a::Int, b::Int)
    length(regions) < 2 && return regions
    best = -1; best2 = -1; best_i = 0
    for (i, r) in enumerate(regions)
        r.extra === nothing && continue
        if r.extra.dp_max > best
            best2 = best; best = r.extra.dp_max; best_i = i
        elseif r.extra.dp_max > best2
            best2 = r.extra.dp_max
        end
    end
    (best_i == 0 || best < 0 || best2 < 0) && return regions
    regions[best_i].qe - regions[best_i].qs < Float64(qlen) * frac && return regions
    best2 < Float64(best) * frac && return regions
    div = 1.0 - event_identity(regions[best_i])
    div < 0.02 && (div = 0.02)
    b2 = 0.5 / div
    b2 * a < b && (b2 = Float64(a) / b)
    for r in regions
        r.extra === nothing && continue
        r.extra.dp_max = max(_recal_max_dp(r, b2, a), Int32(0))
    end
    return regions
end

# mm_squeeze_a: keep only the anchors referenced by regions, in region order of `as`.
function _squeeze_anchors!(regions::Vector{Region}, a::Vector{MM128})
    order = sortperm([Int(r.as) for r in regions])
    out = MM128[]
    for i in order
        r = regions[i]
        start = Int(r.as)
        r.as = length(out)
        append!(out, @view a[(start + 1):(start + r.cnt)])
    end
    return out
end

"""
    align_regions!(regions, index, sequence, a, opt) -> regions

Base-level alignment of the chained regions (minimap2 `mm_align_skeleton`), then
filtering and sorting by DP score.
"""
function align_regions!(regions::Vector{Region}, index::MinimapIndex, sequence::AbstractVector{UInt8},
                        a::Vector{MM128}, opt::MapOptions)
    qlen = length(sequence)
    forward = UInt8[nt4(b) for b in sequence]
    reverse_codes = UInt8[c < 4 ? 0x03 - c : 0x04 for c in Iterators.reverse(forward)]
    qseqs = (forward, reverse_codes)
    a = _squeeze_anchors!(regions, a)
    ez = KswResult()
    i = 1
    while i <= length(regions)
        r2 = _align1!(opt, index, qlen, qseqs, regions[i], a, ez)
        r2 === nothing || insert!(regions, i + 1, r2)
        if i > 1 && regions[i].split_inv && (opt.flag & F_NO_INV) == 0
            r_inv = _align1_inv(opt, index, qlen, qseqs, regions[i - 1], regions[i], ez)
            if r_inv !== nothing
                insert!(regions, i + 1, r_inv)
                i += 1
            end
        end
        i += 1
    end
    filter_regions!(regions, opt, qlen)
    if (opt.flag & (F_SR | F_SR_RNA | F_ALL_CHAINS)) == 0 && qlen >= opt.rank_min_len
        _update_dp_max!(qlen, regions, opt.rank_frac, opt.a, opt.b)
        filter_regions!(regions, opt, qlen)
    end
    hit_sort!(regions, opt.alt_drop)
    return regions, a
end
