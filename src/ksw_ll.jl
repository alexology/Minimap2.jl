# Striped Smith–Waterman score (ksw2_ll_sse.c), used by minimap2 to test for
# inversions.

"""
    ksw_ll_i16(query, target, m, mat, gapo, gape) -> (score, qe, te)

Local alignment score with 0-based end positions on the query and the target
(ksw2 `ksw_ll_i16`). Upstream's striped Smith–Waterman (8 lanes of 16 bits with a
lazy F loop) is emulated literally: as upstream, the column maximum ignores
lazy-F corrections, the query is padded to a multiple of 8 with zero scores, and
`qe` is the last position, in striped memory order, reaching the best score in
the last target column that reached it.
"""
function ksw_ll_i16(query::AbstractVector{UInt8}, target::AbstractVector{UInt8}, m::Int, mat::Vector{Int8},
                    gapo::Int, gape::Int)
    qlen = length(query)
    slen = (qlen + 7) ÷ 8
    # Striped vectors: segment j, lane l holds query position l * slen + j.
    profile = zeros(Int, 8, slen, m)
    for a in 0:(m - 1), j in 0:(slen - 1), l in 0:7
        k = j + l * slen
        profile[l + 1, j + 1, a + 1] = k >= qlen ? 0 : Int(mat[a * m + Int(query[k + 1]) + 1])
    end
    H0 = zeros(Int, 8, slen); H1 = zeros(Int, 8, slen); E = zeros(Int, 8, slen); Hmax = zeros(Int, 8, slen)
    gapoe = gapo + gape
    sat16(x) = clamp(x, -32768, 32767)
    subs_u16(x, y) = max(x - y, 0)          # operands are non-negative here
    gmax = 0; te = -1
    h = zeros(Int, 8); f = zeros(Int, 8); vmax = zeros(Int, 8)
    for i in 0:(length(target) - 1)
        a = Int(target[i + 1])
        fill!(f, 0); fill!(vmax, 0)
        h[1] = 0
        for l in 2:8
            h[l] = H0[l - 1, slen]
        end
        for j in 1:slen
            for l in 1:8
                hv = sat16(h[l] + profile[l, j, a + 1])
                ev = E[l, j]
                hv = max(hv, ev); hv = max(hv, f[l])
                vmax[l] = max(vmax[l], hv)
                H1[l, j] = hv
                hv = subs_u16(hv, gapoe)
                ev = max(subs_u16(ev, gape), hv)
                E[l, j] = ev
                f[l] = max(subs_u16(f[l], gape), hv)
                h[l] = H0[l, j]
            end
        end
        done = false
        for _ in 1:8
            for l in 8:-1:2
                f[l] = f[l - 1]
            end
            f[1] = 0
            for j in 1:slen
                any_greater = false
                for l in 1:8
                    hv = max(H1[l, j], f[l])
                    H1[l, j] = hv
                    hv = subs_u16(hv, gapoe)
                    f[l] = subs_u16(f[l], gape)
                    f[l] > hv && (any_greater = true)
                end
                if !any_greater
                    done = true
                    break
                end
            end
            done && break
        end
        imax = maximum(vmax)
        if imax >= gmax
            gmax = imax; te = i
            copyto!(Hmax, H1)
        end
        H0, H1 = H1, H0
    end
    qe = -1
    for i in 0:(slen * 8 - 1)
        j = i ÷ 8; l = i % 8
        Hmax[l + 1, j + 1] == gmax && (qe = j + l * slen)
    end
    return gmax, qe, te
end
