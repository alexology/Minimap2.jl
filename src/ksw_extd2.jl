# ksw2 `ksw_extd2_sse`: extension alignment with dual affine gaps, computed by
# anti-diagonals with 8-bit difference scores. Each cell is computed exactly as one
# lane of the SSE code (wrapping Int8 arithmetic, same tie-breaks). The row loop
# reads the previous anti-diagonal shifted by one cell from copies made before the
# row is updated, which removes the carry between cells and lets it vectorise.
#
# Upstream quirks reproduced: cells are computed in whole 16-cell blocks beyond the
# band, the score buffer keeps scores from earlier rows for those cells, 16-byte
# score stores may run past the score buffer into the target buffer, and loads of
# the target may run into the reversed query that follows it in memory.

"""
    ksw_extd2!(ez, query, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag) -> ez

Extension alignment with dual affine gaps (ksw2 `ksw_extd2_sse`); `query` and
`target` hold codes `0:m-1`, `mat` is the `m × m` score matrix in row-major order.
"""
function ksw_extd2!(ez::KswResult, query::AbstractVector{UInt8}, target::AbstractVector{UInt8}, m::Int,
                    mat::Vector{Int8}, q::Int, e::Int, q2::Int, e2::Int, w::Int, zdrop::Int, end_bonus::Int, flag::Int)
    _reset!(ez)
    (m <= 1 || isempty(query) || isempty(target)) && return ez
    # Make sure q + e is no larger than q2 + e2 (the core must not reassign q and e:
    # they are captured by closures).
    q2 + e2 < q + e && return _ksw_extd2_core!(ez, query, target, m, mat, q2, e2, q, e, w, zdrop, end_bonus, flag)
    return _ksw_extd2_core!(ez, query, target, m, mat, q, e, q2, e2, w, zdrop, end_bonus, flag)
end

function _ksw_extd2_core!(ez::KswResult, query::AbstractVector{UInt8}, target::AbstractVector{UInt8}, m::Int,
                          mat::Vector{Int8}, q::Int, e::Int, q2::Int, e2::Int, w::Int, zdrop::Int, end_bonus::Int, flag::Int)
    qlen = length(query); tlen = length(target)
    with_cigar = (flag & KSW_EZ_SCORE_ONLY) == 0
    approx_max = (flag & KSW_EZ_APPROX_MAX) != 0
    right = (flag & KSW_EZ_RIGHT) != 0
    qe = q + e
    sc_mch = mat[1]; sc_mis = mat[2]
    sc_N = mat[m * m] == 0 ? Int8(-e2) : mat[m * m]
    m1 = UInt8(m - 1)
    wl = wr = w < 0 ? max(tlen, qlen) : w
    tlen_ = (tlen + 15) ÷ 16
    n_col_ = min(qlen, tlen)
    n_col_ = (min(n_col_, wl + 1) + 15) ÷ 16 + 1
    qlen_ = (qlen + 15) ÷ 16
    min_sc = Int(mat[2])
    for t in 2:(m * m)
        min_sc = min(min_sc, Int(mat[t]))
    end
    -min_sc > 2 * (q + e) && return ez
    long_thres = e != e2 ? (q2 - q) ÷ (e - e2) - 1 : 0
    (q2 + e2 + long_thres * e2 > q + e + long_thres * e) && (long_thres += 1)
    long_diff = long_thres * (e - e2) - (q2 - q) - e2
    seg = tlen_ * 16; cap_qr = (qlen_ + 1) * 16
    u, v, x, y, x2, y2, s, xp, vp, x2p = ez.work8
    foreach(b -> resize!(b, seg), ez.work8)
    fill!(u, Int8(-q - e)); fill!(v, Int8(-q - e)); fill!(x, Int8(-q - e)); fill!(y, Int8(-q - e))
    fill!(x2, Int8(-q2 - e2)); fill!(y2, Int8(-q2 - e2)); fill!(s, Int8(0))
    sf = resize!(ez.sf, seg); fill!(sf, 0x00); copyto!(sf, 1, target, 1, tlen)
    qr = resize!(ez.qr, cap_qr); fill!(qr, 0x00)
    for t in 0:(qlen - 1)
        qr[t + 1] = query[qlen - t]
    end
    H = ez.H
    approx_max || (resize!(H, seg); fill!(H, Int32(KSW_NEG_INF)))
    nrow = qlen + tlen - 1
    n_col = n_col_ * 16
    # Every cell the backtrack reads is written first, so p needs no clearing.
    p = with_cigar ? resize!(ez.p, n_col * nrow) : ez.p
    off = resize!(ez.off, with_cigar ? nrow : 0); off_end = resize!(ez.off_end, with_cigar ? nrow : 0)
    H0 = Int32(0); last_H0_t = 0
    last_st = last_en = -1
    q8 = Int8(q); q28 = Int8(q2); qe8 = Int8(qe); qe28 = Int8(q2 + e2); zero8 = Int8(0)
    # Byte loads/stores past the end of a buffer continue in the next buffer, as in
    # upstream's single allocation (s, then sf, then qr).
    sf_load(i) = i < seg ? sf[i + 1] : (i - seg < cap_qr ? qr[i - seg + 1] : 0x00)
    qr_load(i) = i < cap_qr ? qr[i + 1] : 0x00
    function s_store!(i, value::Int8)
        if i < seg
            s[i + 1] = value
        elseif i - seg < seg
            sf[i - seg + 1] = reinterpret(UInt8, value)
        end
        return nothing
    end
    chunk = Vector{Int8}(undef, 16)
    for r in 0:(nrow - 1)
        st = 0; en = tlen - 1
        st < r - qlen + 1 && (st = r - qlen + 1)
        en > r && (en = r)
        st < (r - wr + 1) >> 1 && (st = (r - wr + 1) >> 1)
        en > (r + wl) >> 1 && (en = (r + wl) >> 1)
        if st > en
            ez.zdropped = true
            break
        end
        st0 = st; en0 = en
        st = st ÷ 16 * 16; en = (en + 16) ÷ 16 * 16 - 1
        if st > 0
            if st - 1 >= last_st && st - 1 <= last_en
                x1 = x[st]; x21 = x2[st]; v1 = v[st]
            else
                x1 = Int8(-q - e); x21 = Int8(-q2 - e2); v1 = Int8(-q - e)
            end
        else
            x1 = Int8(-q - e); x21 = Int8(-q2 - e2)
            v1 = Int8(r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2)
        end
        if en >= r
            y[r + 1] = Int8(-q - e); y2[r + 1] = Int8(-q2 - e2)
            u[r + 1] = Int8(r == 0 ? -q - e : r < long_thres ? -e : r == long_thres ? long_diff : -e2)
        end
        qrr = qlen - 1 - r
        last_chunk = st0 + ((en0 - st0) ÷ 16 + 1) * 16 - 1    # last byte written by the 16-byte stores
        if (flag & KSW_EZ_GENERIC_SC) == 0 && last_chunk < seg && qrr + last_chunk < cap_qr
            # The whole row lies inside the buffers: the chunked stores reduce to one loop.
            @inbounds @simd for i in (st0 + 1):(last_chunk + 1)
                sq = sf[i]; sb = qr[qrr + i]
                value = ifelse(sq == sb, sc_mch, sc_mis)
                s[i] = ifelse((sq == m1) | (sb == m1), sc_N, value)
            end
        elseif (flag & KSW_EZ_GENERIC_SC) == 0
            t = st0
            @inbounds while t <= en0
                if t + 15 < seg && qrr + t + 15 < cap_qr
                    for k in 1:16
                        sq = sf[t + k]; sb = qr[qrr + t + k]
                        value = ifelse(sq == sb, sc_mch, sc_mis)
                        chunk[k] = ifelse((sq == m1) | (sb == m1), sc_N, value)
                    end
                    copyto!(s, t + 1, chunk, 1, 16)
                else
                    for k in 0:15
                        sq = sf_load(t + k); sb = qr_load(qrr + t + k)
                        value = sq == sb ? sc_mch : sc_mis
                        (sq == m1 || sb == m1) && (value = sc_N)
                        chunk[k + 1] = value
                    end
                    for k in 0:15
                        s_store!(t + k, chunk[k + 1])
                    end
                end
                t += 16
            end
        else
            for t in st0:en0
                s_store!(t, mat[Int(sf_load(t)) * m + Int(qr_load(qrr + t)) + 1])
            end
        end
        # Previous anti-diagonal shifted by one cell: xp[t] = x[t - 1], with the
        # boundary values at t = st.
        n = en - st + 1
        xp[st + 1] = x1; vp[st + 1] = v1; x2p[st + 1] = x21
        if n > 1
            copyto!(xp, st + 2, x, st + 1, n - 1)
            copyto!(vp, st + 2, v, st + 1, n - 1)
            copyto!(x2p, st + 2, x2, st + 1, n - 1)
        end
        if !with_cigar
            @inbounds @simd for i in (st + 1):(en + 1)
                z = s[i]; vt1 = vp[i]; ut = u[i]
                a = xp[i] + vt1; b = y[i] + ut; a2 = x2p[i] + vt1; b2 = y2[i] + ut
                z = min(max(max(max(max(z, a), b), a2), b2), sc_mch)
                u[i] = z - vt1; v[i] = z - ut
                tmp = z - q8; a -= tmp; b -= tmp
                tmp = z - q28; a2 -= tmp; b2 -= tmp
                x[i] = max(a, zero8) - qe8; y[i] = max(b, zero8) - qe8
                x2[i] = max(a2, zero8) - qe28; y2[i] = max(b2, zero8) - qe28
            end
        else
            off[r + 1] = st; off_end[r + 1] = en
            column = r * n_col - st
            if !right
                @inbounds @simd for i in (st + 1):(en + 1)
                    z = s[i]; vt1 = vp[i]; ut = u[i]
                    a = xp[i] + vt1; b = y[i] + ut; a2 = x2p[i] + vt1; b2 = y2[i] + ut
                    d = ifelse(a > z, 0x01, 0x00); z = max(z, a)
                    d = ifelse(b > z, 0x02, d); z = max(z, b)
                    d = ifelse(a2 > z, 0x03, d); z = max(z, a2)
                    d = ifelse(b2 > z, 0x04, d); z = max(z, b2)
                    z = min(z, sc_mch)
                    u[i] = z - vt1; v[i] = z - ut
                    tmp = z - q8; a -= tmp; b -= tmp
                    tmp = z - q28; a2 -= tmp; b2 -= tmp
                    x[i] = max(a, zero8) - qe8; y[i] = max(b, zero8) - qe8
                    x2[i] = max(a2, zero8) - qe28; y2[i] = max(b2, zero8) - qe28
                    d |= ifelse(a > zero8, 0x08, 0x00) | ifelse(b > zero8, 0x10, 0x00) |
                         ifelse(a2 > zero8, 0x20, 0x00) | ifelse(b2 > zero8, 0x40, 0x00)
                    p[column + i] = d
                end
            else
                @inbounds @simd for i in (st + 1):(en + 1)
                    z = s[i]; vt1 = vp[i]; ut = u[i]
                    a = xp[i] + vt1; b = y[i] + ut; a2 = x2p[i] + vt1; b2 = y2[i] + ut
                    d = ifelse(z > a, 0x00, 0x01); z = max(z, a)
                    d = ifelse(z > b, d, 0x02); z = max(z, b)
                    d = ifelse(z > a2, d, 0x03); z = max(z, a2)
                    d = ifelse(z > b2, d, 0x04); z = max(z, b2)
                    z = min(z, sc_mch)
                    u[i] = z - vt1; v[i] = z - ut
                    tmp = z - q8; a -= tmp; b -= tmp
                    tmp = z - q28; a2 -= tmp; b2 -= tmp
                    x[i] = max(a, zero8) - qe8; y[i] = max(b, zero8) - qe8
                    x2[i] = max(a2, zero8) - qe28; y2[i] = max(b2, zero8) - qe28
                    d |= ifelse(a >= zero8, 0x08, 0x00) | ifelse(b >= zero8, 0x10, 0x00) |
                         ifelse(a2 >= zero8, 0x20, 0x00) | ifelse(b2 >= zero8, 0x40, 0x00)
                    p[column + i] = d
                end
            end
        end
        if !approx_max
            @inbounds if r > 0
                en1 = st0 + (en0 - st0) ÷ 4 * 4
                H[en0 + 1] = en0 > 0 ? H[en0] + Int32(u[en0 + 1]) : H[en0 + 1] + Int32(v[en0 + 1])
                max_H = H[en0 + 1]; max_t = en0
                # SSE: four lanes over [st0, en1), lane k holding positions st0 + k + 4j, each
                # keeping its first strictly larger value, then reduced in lane order. The
                # winner is the first occurrence of the overall maximum in the lowest lane.
                @simd for i in (st0 + 1):en1
                    H[i] += Int32(v[i])
                end
                if en1 > st0
                    best = maximum(view(H, (st0 + 1):en1))
                    if best > max_H
                        best_lane = 4; best_t = 0
                        for i in (st0 + 1):en1
                            if H[i] == best
                                lane = (i - 1 - st0) & 3
                                if lane < best_lane
                                    best_lane = lane; best_t = i - 1
                                    lane == 0 && break
                                end
                            end
                        end
                        max_H = best; max_t = best_t
                    end
                end
                t = en1
                while t < en0
                    H[t + 1] += Int32(v[t + 1])
                    if H[t + 1] > max_H
                        max_H = H[t + 1]; max_t = t
                    end
                    t += 1
                end
            else
                H[1] = Int32(v[1]) - Int32(qe); max_H = H[1]; max_t = 0
            end
            if en0 == tlen - 1 && H[en0 + 1] > ez.mte
                ez.mte = H[en0 + 1]; ez.mte_q = r - en0
            end
            if r - st0 == qlen - 1 && H[st0 + 1] > ez.mqe
                ez.mqe = H[st0 + 1]; ez.mqe_t = st0
            end
            _apply_zdrop!(ez, max_H, r, max_t, zdrop, e2) && break
            (r == qlen + tlen - 2 && en0 == tlen - 1) && (ez.score = H[tlen])
        else
            if r > 0
                if last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0
                    d0 = Int32(v[last_H0_t + 1]); d1 = Int32(u[last_H0_t + 2])
                    if d0 > d1
                        H0 += d0
                    else
                        H0 += d1; last_H0_t += 1
                    end
                elseif last_H0_t >= st0 && last_H0_t <= en0
                    H0 += Int32(v[last_H0_t + 1])
                else
                    last_H0_t += 1; H0 += Int32(u[last_H0_t + 1])
                end
            else
                H0 = Int32(v[1]) - Int32(qe); last_H0_t = 0
            end
            ((flag & KSW_EZ_APPROX_DROP) != 0 && _apply_zdrop!(ez, H0, r, last_H0_t, zdrop, e2)) && break
            (r == qlen + tlen - 2 && en0 == tlen - 1) && (ez.score = H0)
        end
        last_st = st; last_en = en
    end
    if with_cigar
        rev_cigar = (flag & KSW_EZ_REV_CIGAR) != 0
        if !ez.zdropped && (flag & KSW_EZ_EXTZ_ONLY) == 0
            ez.cigar = _ksw_backtrack(rev_cigar, p, n_col, off, off_end, tlen - 1, qlen - 1)
        elseif !ez.zdropped && (flag & KSW_EZ_EXTZ_ONLY) != 0 && ez.mqe + end_bonus > ez.max
            ez.reach_end = true
            ez.cigar = _ksw_backtrack(rev_cigar, p, n_col, off, off_end, Int(ez.mqe_t), qlen - 1)
        elseif ez.max_t >= 0 && ez.max_q >= 0
            ez.cigar = _ksw_backtrack(rev_cigar, p, n_col, off, off_end, Int(ez.max_t), Int(ez.max_q))
        end
    end
    return ez
end
