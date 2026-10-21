# minimap2's sorting primitives (ksort.h, misc.c). The in-place radix sort is not
# stable, so the order of equal keys is reproduced by porting it literally.

"""
    MM128

minimap2's `mm128_t`: a pair of 64-bit words `x` and `y` whose meaning depends on
the context (minimizers, anchors, ...).
"""
struct MM128
    x::UInt64
    y::UInt64
end

const _RS_MIN_SIZE = 64
const _RS_MAX_BITS = 8

function _insertion_sort!(v::Vector{T}, lo::Int, hi::Int, key) where {T}
    # Sorts v[lo:hi] (rs_insertsort).
    for i in (lo + 1):hi
        if key(v[i]) < key(v[i - 1])
            tmp = v[i]
            j = i
            while j > lo && key(tmp) < key(v[j - 1])
                v[j] = v[j - 1]
                j -= 1
            end
            v[j] = tmp
        end
    end
    return v
end

function _radix_pass!(v::Vector{T}, lo::Int, hi::Int, n_bits::Int, s::Int, key) where {T}
    # rs_sort over v[lo:hi]; buckets hold 1-based half-open [b, e).
    size = 1 << n_bits; m = UInt64(size - 1)
    b = fill(lo, size); e = fill(lo, size)
    for i in lo:hi
        e[Int((key(v[i]) >> s) & m) + 1] += 1
    end
    for k in 2:size
        e[k] += e[k - 1] - lo
        b[k] = e[k - 1]
    end
    k = 1
    while k <= size
        if b[k] != e[k]
            l = Int((key(v[b[k]]) >> s) & m) + 1
            if l != k
                tmp = v[b[k]]
                while true
                    swap = tmp; tmp = v[b[l]]; v[b[l]] = swap; b[l] += 1
                    l = Int((key(tmp) >> s) & m) + 1
                    l == k && break
                end
                v[b[k]] = tmp; b[k] += 1
            else
                b[k] += 1
            end
        else
            k += 1
        end
    end
    b[1] = lo
    for k in 2:size
        b[k] = e[k - 1]
    end
    if s != 0
        s = s > n_bits ? s - n_bits : 0
        for k in 1:size
            if e[k] - b[k] > _RS_MIN_SIZE
                _radix_pass!(v, b[k], e[k] - 1, n_bits, s, key)
            elseif e[k] - b[k] > 1
                _insertion_sort!(v, b[k], e[k] - 1, key)
            end
        end
    end
    return v
end

# radix_sort_128x / radix_sort_64 on v[lo:hi].
function _radix_sort!(v::Vector{T}, lo::Int=1, hi::Int=length(v); key=identity) where {T}
    if hi - lo + 1 <= _RS_MIN_SIZE
        _insertion_sort!(v, lo, hi, key)
    else
        _radix_pass!(v, lo, hi, _RS_MAX_BITS, 7 * _RS_MAX_BITS, key)
    end
    return v
end
radix_sort_128x!(v::Vector{MM128}, lo::Int=1, hi::Int=length(v)) = _radix_sort!(v, lo, hi; key=a -> a.x)

# ksort heap on a 1-based vector of length n (ks_heapdown, ks_heapmake), max-heap
# under `lt`.
function _heapdown!(l::AbstractVector, i::Int, n::Int, lt)
    # i is 0-based as in ksort.
    k = i
    tmp = l[i + 1]
    while (k = (k << 1) + 1) < n
        if k != n - 1 && lt(l[k + 1], l[k + 2])
            k += 1
        end
        lt(l[k + 1], tmp) && break
        l[i + 1] = l[k + 1]; i = k
    end
    l[i + 1] = tmp
    return l
end

function _heapmake!(l::AbstractVector, n::Int, lt)
    for i in ((n >> 1) - 1):-1:0
        _heapdown!(l, i, n, lt)
    end
    return l
end
