# Chaining (lchain.c): dynamic programming and RMQ-based chaining, backtracking.
# Chains are returned as `u` entries (score << 32 | anchor count) with the anchors
# of each chain stored consecutively, sorted by reference position.

# mg_log2 (mmpriv.h): fast float log2, valid for x >= 2.
@inline function _mg_log2(x::Float32)
    bits = reinterpret(UInt32, x)
    log2_value = Float32(Int32((bits >> 23) & 0xff) - Int32(128))
    bits &= ~(UInt32(255) << 23)
    bits += UInt32(127) << 23
    z = reinterpret(Float32, bits)
    log2_value += (-0.34484843f0 * z + 2.02466578f0) * z - 0.67487759f0
    return log2_value
end

@inline _ay(a::MM128) = reinterpret(Int32, UInt32(a.y & 0xffffffff))
@inline _ax(a::MM128) = reinterpret(Int32, UInt32(a.x & 0xffffffff))
@inline _span(a::MM128) = Int32((a.y >> 32) & 0xff)
@inline _seg(a::MM128) = Int32((a.y & SEED_SEG_MASK) >> SEED_SEG_SHIFT)

function _chain_score(ai::MM128, aj::MM128, max_dist_x::Int32, max_dist_y::Int32, bw::Int32,
                      pen_gap::Float32, pen_skip::Float32, is_cdna::Bool, n_seg::Int)
    @inbounds begin
    dq = _ay(ai) - _ay(aj)
    (dq <= 0 || dq > max_dist_x) && return typemin(Int32)
    dr = reinterpret(Int32, UInt32((ai.x - aj.x) & 0xffffffff))
    sidi = _seg(ai); sidj = _seg(aj)
    (sidi == sidj && (dr == 0 || dq > max_dist_y)) && return typemin(Int32)
    dd = dr > dq ? dr - dq : dq - dr
    (sidi == sidj && dd > bw) && return typemin(Int32)
    (n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) && return typemin(Int32)
    dg = dr < dq ? dr : dq
    q_span = _span(aj)
    sc = q_span < dg ? q_span : dg
    if dd != 0 || dg > q_span
        lin_pen = pen_gap * Float32(dd) + pen_skip * Float32(dg)
        log_pen = dd >= 1 ? _mg_log2(Float32(dd + 1)) : 0.0f0
        if is_cdna || sidi != sidj
            if sidi != sidj && dr == 0
                sc += Int32(1)
            elseif dr > dq || sidi != sidj
                sc -= unsafe_trunc(Int32, lin_pen < log_pen ? lin_pen : log_pen)
            else
                sc -= unsafe_trunc(Int32, lin_pen + 0.5f0 * log_pen)
            end
        else
            sc -= unsafe_trunc(Int32, lin_pen + 0.5f0 * log_pen)
        end
    end
    return sc
    end
end

# mg_chain_bk_end; z holds (score, index) pairs, 0-based indices, p[i] = -1 for none.
function _chain_bk_end(max_drop::Int32, z::Vector{MM128}, f::Vector{Int32}, p::Vector{Int64}, t::Vector{Int32}, k::Int)
    @inbounds begin
    i = Int64(z[k].y); end_i = Int64(-1); max_i = i
    max_s = Int32(0)
    (i < 0 || t[i + 1] != 0) && return i
    while true
        t[i + 1] = 2
        end_i = i = p[i + 1]
        s = i < 0 ? reinterpret(Int32, UInt32(z[k].x & 0xffffffff)) : reinterpret(Int32, UInt32(z[k].x & 0xffffffff)) - f[i + 1]
        if s > max_s
            max_s = s; max_i = i
        elseif max_s - s > max_drop
            break
        end
        (i >= 0 && t[i + 1] == 0) || break
    end
    i = Int64(z[k].y)
    while i >= 0 && i != end_i
        t[i + 1] = 0
        i = p[i + 1]
    end
    return max_i
    end
end

# mg_chain_backtrack: returns (u, v) with v the anchor indices of each chain from
# its last anchor backwards.
function _chain_backtrack(n::Int, f::Vector{Int32}, p::Vector{Int64}, t::Vector{Int32},
                          min_cnt::Int, min_sc::Int, max_drop::Int32)
    @inbounds begin
    z = MM128[]
    for i in 0:(n - 1)
        f[i + 1] >= min_sc && push!(z, MM128(reinterpret(UInt32, f[i + 1]) % UInt64, UInt64(i)))
    end
    isempty(z) && return UInt64[], Int32[]
    radix_sort_128x!(z)
    u = UInt64[]; v = Int32[]
    fill!(t, 0)
    for k in length(z):-1:1
        t[z[k].y + 1] == 0 || continue
        n_v0 = length(v)
        end_i = _chain_bk_end(max_drop, z, f, p, t, k)
        i = Int64(z[k].y)
        while i != end_i
            push!(v, Int32(i)); t[i + 1] = 1
            i = p[i + 1]
        end
        sc = i < 0 ? reinterpret(Int32, UInt32(z[k].x & 0xffffffff)) : reinterpret(Int32, UInt32(z[k].x & 0xffffffff)) - f[i + 1]
        if sc >= min_sc && length(v) > n_v0 && length(v) - n_v0 >= min_cnt
            push!(u, UInt64(reinterpret(UInt32, sc)) << 32 | UInt64(length(v) - n_v0))
        else
            resize!(v, n_v0)
        end
    end
    return u, v
    end
end

# compact_a: anchors of each chain in forward order, chains sorted by the
# reference position of their first anchor.
function _compact_anchors(u::Vector{UInt64}, v::Vector{Int32}, a::Vector{MM128})
    b = MM128[]
    k = 0
    for ui in u
        ni = Int(ui & 0xffffffff)
        for j in 0:(ni - 1)
            push!(b, a[v[k + (ni - j - 1) + 1] + 1])
        end
        k += ni
    end
    w = Vector{MM128}(undef, length(u))
    k = 0
    for i in eachindex(u)
        w[i] = MM128(b[k + 1].x, UInt64(k) << 32 | UInt64(i - 1))
        k += Int(u[i] & 0xffffffff)
    end
    radix_sort_128x!(w)
    u2 = similar(u); out = MM128[]
    for i in eachindex(w)
        j = Int(w[i].y & 0xffffffff) + 1; n = Int(u[j] & 0xffffffff)
        u2[i] = u[j]
        start = Int(w[i].y >> 32)
        append!(out, @view b[(start + 1):(start + n)])
    end
    return u2, out
end

"""
    chain_dp(anchors; max_dist_x, max_dist_y, bw, max_skip, max_iter, min_cnt, min_sc,
             pen_gap, pen_skip, is_cdna=false, n_seg=1) -> (u, anchors)

minimap2's DP chaining (`mg_lchain_dp`).
"""
function chain_dp(a::Vector{MM128}, max_dist_x::Int, max_dist_y::Int, bw::Int, max_skip::Int, max_iter::Int,
                  min_cnt::Int, min_sc::Int, pen_gap::Float32, pen_skip::Float32, is_cdna::Bool, n_seg::Int)
    @inbounds begin
    n = length(a)
    n == 0 && return UInt64[], MM128[]
    max_dist_x < bw && (max_dist_x = bw)
    (max_dist_y < bw && !is_cdna) && (max_dist_y = bw)
    max_drop = is_cdna ? typemax(Int32) : Int32(bw)
    mdx = Int32(max_dist_x); mdy = Int32(max_dist_y); bw32 = Int32(bw)
    p = Vector{Int64}(undef, n); f = Vector{Int32}(undef, n)
    v = Vector{Int32}(undef, n); t = zeros(Int32, n)
    st = 0; max_ii = -1
    for i in 0:(n - 1)
        ai = a[i + 1]
        max_j = -1
        max_f = _span(ai); n_skip = 0
        while st < i && ((ai.x >> 32) != (a[st + 1].x >> 32) || ai.x > a[st + 1].x + UInt64(max_dist_x))
            st += 1
        end
        i - st > max_iter && (st = i - max_iter)
        j = i - 1
        while j >= st
            sc = _chain_score(ai, a[j + 1], mdx, mdy, bw32, pen_gap, pen_skip, is_cdna, n_seg)
            if sc != typemin(Int32)
                sc += f[j + 1]
                if sc > max_f
                    max_f = sc; max_j = j
                    n_skip > 0 && (n_skip -= 1)
                elseif t[j + 1] == Int32(i)
                    n_skip += 1
                    n_skip > max_skip && break
                end
                p[j + 1] >= 0 && (t[p[j + 1] + 1] = Int32(i))
            end
            j -= 1
        end
        end_j = j
        if max_ii < 0 || ai.x - a[max_ii + 1].x > UInt64(max_dist_x)   # unsigned, as in C
            best = typemin(Int32); max_ii = -1
            for jj in (i - 1):-1:st
                if best < f[jj + 1]
                    best = f[jj + 1]; max_ii = jj
                end
            end
        end
        if max_ii >= 0 && max_ii < end_j
            tmp = _chain_score(ai, a[max_ii + 1], mdx, mdy, bw32, pen_gap, pen_skip, is_cdna, n_seg)
            if tmp != typemin(Int32) && max_f < tmp + f[max_ii + 1]
                max_f = tmp + f[max_ii + 1]; max_j = max_ii
            end
        end
        f[i + 1] = max_f; p[i + 1] = max_j
        v[i + 1] = max_j >= 0 && v[max_j + 1] > max_f ? v[max_j + 1] : max_f
        if max_ii < 0 || (ai.x - a[max_ii + 1].x <= UInt64(max_dist_x) && f[max_ii + 1] < f[i + 1])
            max_ii = i
        end
    end
    u, vv = _chain_backtrack(n, f, p, t, min_cnt, min_sc, max_drop)
    isempty(u) && return UInt64[], MM128[]
    return _compact_anchors(u, vv, a)
    end
end

function _chain_score_simple(ai::MM128, aj::MM128, pen_gap::Float32, pen_skip::Float32)
    @inbounds begin
    dq = _ay(ai) - _ay(aj)
    dr = reinterpret(Int32, UInt32((ai.x - aj.x) & 0xffffffff))
    dd = dr > dq ? dr - dq : dq - dr
    width = dd
    dg = dr < dq ? dr : dq
    q_span = _span(aj)
    sc = q_span < dg ? q_span : dg
    exact = dd == 0 && dg <= q_span
    if dd != 0 || dq > q_span
        lin_pen = pen_gap * Float32(dd) + pen_skip * Float32(dg)
        log_pen = dd >= 1 ? _mg_log2(Float32(dd + 1)) : 0.0f0
        sc -= unsafe_trunc(Int32, lin_pen + 0.5f0 * log_pen)
    end
    return sc, exact, width
    end
end

"""
    chain_rmq(anchors, max_dist, max_dist_inner, bw, max_chn_skip, cap_rmq_size, min_cnt,
              min_sc, pen_gap, pen_skip) -> (u, anchors)

minimap2's RMQ-based chaining (`mg_lchain_rmq`).
"""
function chain_rmq(a::Vector{MM128}, max_dist::Int, max_dist_inner::Int, bw::Int, max_chn_skip::Int,
                   cap_rmq_size::Int, min_cnt::Int, min_sc::Int, pen_gap::Float32, pen_skip::Float32)
    @inbounds begin
    n = length(a)
    n == 0 && return UInt64[], MM128[]
    max_dist < bw && (max_dist = bw)
    max_dist_inner < 0 && (max_dist_inner = 0)
    max_dist_inner > max_dist && (max_dist_inner = max_dist)
    max_drop = Int32(bw)
    p = Vector{Int64}(undef, n); f = Vector{Int32}(undef, n)
    t = zeros(Int32, n); v = Vector{Int32}(undef, n)
    pool = _task_rmq_pool()
    root = Ref(0); root_inner = Ref(0)
    st = 0; st_inner = 0; i0 = 0
    for i in 0:(n - 1)
        ai = a[i + 1]
        max_j = -1
        q_span = _span(ai); max_f = q_span
        if i0 < i && a[i0 + 1].x != ai.x
            for j in i0:(i - 1)
                aj = a[j + 1]
                pri = -(Float64(f[j + 1]) + 0.5 * Float64(pen_gap) * Float64(Int64(_ax(aj)) + Int64(_ay(aj))))
                q = _rmq_alloc!(pool, _ay(aj), Int64(j), pri)
                _rmq_insert!(pool, root, q)
                if max_dist_inner > 0
                    r = _rmq_alloc!(pool, _ay(aj), Int64(j), pri)
                    _rmq_insert!(pool, root_inner, r)
                end
            end
            i0 = i
        end
        while st < i && ((ai.x >> 32) != (a[st + 1].x >> 32) || ai.x > a[st + 1].x + UInt64(max_dist) ||
                         _rmq_size(pool, root[]) > cap_rmq_size)
            q = _rmq_find(pool, root[], _ay(a[st + 1]), Int64(st))
            q != 0 && _rmq_erase!(pool, root, q)
            st += 1
        end
        if max_dist_inner > 0
            while st_inner < i && ((ai.x >> 32) != (a[st_inner + 1].x >> 32) || ai.x > a[st_inner + 1].x + UInt64(max_dist_inner) ||
                                   _rmq_size(pool, root_inner[]) > cap_rmq_size)
                q = _rmq_find(pool, root_inner[], _ay(a[st_inner + 1]), Int64(st_inner))
                q != 0 && _rmq_erase!(pool, root_inner, q)
                st_inner += 1
            end
        end
        q = _rmq_query(pool, root[], _ay(ai) - Int32(max_dist), Int64(typemax(Int32)), _ay(ai), Int64(0))
        if q != 0
            n_skip = 0
            j = Int(pool.i[q])
            sc0, exact, width = _chain_score_simple(ai, a[j + 1], pen_gap, pen_skip)
            sc = f[j + 1] + sc0
            if width <= bw && sc > max_f
                max_f = sc; max_j = j
            end
            if !exact && root_inner[] != 0 && _ay(ai) > 0
                lo, _ = _rmq_interval(pool, root_inner[], _ay(ai) - Int32(1), Int64(n))
                if lo != 0
                    stack = _rmq_itr_find(pool, root_inner[], lo)
                    while !isempty(stack)
                        node = stack[end]
                        pool.y[node] < _ay(ai) - Int32(max_dist_inner) && break
                        j = Int(pool.i[node])
                        sc0, _, width = _chain_score_simple(ai, a[j + 1], pen_gap, pen_skip)
                        sc = f[j + 1] + sc0
                        if width <= bw
                            if sc > max_f
                                max_f = sc; max_j = j
                                n_skip > 0 && (n_skip -= 1)
                            elseif t[j + 1] == Int32(i)
                                n_skip += 1
                                n_skip > max_chn_skip && break
                            end
                            p[j + 1] >= 0 && (t[p[j + 1] + 1] = Int32(i))
                        end
                        _rmq_itr_prev!(pool, stack) || break
                    end
                end
            end
        end
        f[i + 1] = max_f; p[i + 1] = max_j
        v[i + 1] = max_j >= 0 && v[max_j + 1] > max_f ? v[max_j + 1] : max_f
    end
    u, vv = _chain_backtrack(n, f, p, t, min_cnt, min_sc, max_drop)
    isempty(u) && return UInt64[], MM128[]
    return _compact_anchors(u, vv, a)
    end
end
