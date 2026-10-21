# AVL tree with range-minimum queries (krmq.h), specialised for the RMQ chaining
# elements of lchain.c: key (y, i), priority `pri` (smaller is better). Ties in
# priority are resolved by the tree shape, so the tree operations are ported
# literally. Nodes are indices into the pool; 0 is the null node.

mutable struct _RMQPool
    y::Vector{Int32}
    i::Vector{Int64}
    pri::Vector{Float64}
    child::Matrix{Int}      # 2 × capacity: left (1) and right (2)
    s::Vector{Int}          # node with the minimum priority in the subtree
    balance::Vector{Int8}
    size::Vector{Int}
    # Scratch buffers (depth <= 64, as KRMQ_MAX_DEPTH).
    path::Vector{Int}
    dir::Vector{Int}
    path2::Vector{Int}
    cmp1::Vector{Int}
    cmp2::Vector{Int}
    stack::Vector{Int}
    stack_insert::Vector{Int}
end
# Node 1 is the fake node used by erase (its left child is the root).
function _RMQPool()
    pool = _RMQPool(Int32[], Int64[], Float64[], zeros(Int, 2, 256), Int[], Int8[], Int[],
                    zeros(Int, 66), zeros(Int, 66), zeros(Int, 66), zeros(Int, 66), zeros(Int, 66), Int[], zeros(Int, 66))
    _rmq_alloc!(pool, Int32(0), Int64(0), 0.0)
    return pool
end

# Empty the pool, keeping its memory (node 1, the fake node, stays).
function _rmq_reset!(pool::_RMQPool)
    @inbounds begin
    resize!(pool.y, 1); resize!(pool.i, 1); resize!(pool.pri, 1)
    resize!(pool.s, 1); resize!(pool.balance, 1); resize!(pool.size, 1)
    return pool
    end
end

# A pool per task, reused by successive chaining calls.
_task_rmq_pool() = _rmq_reset!(get!(_RMQPool, task_local_storage(), :minimap2_rmq_pool)::_RMQPool)

function _rmq_alloc!(pool::_RMQPool, y::Int32, i::Int64, pri::Float64)
    @inbounds begin
    push!(pool.y, y); push!(pool.i, i); push!(pool.pri, pri)
    n = length(pool.y)
    if size(pool.child, 2) < n
        pool.child = hcat(pool.child, zeros(Int, 2, size(pool.child, 2)))
    end
    pool.child[1, n] = 0; pool.child[2, n] = 0
    push!(pool.s, n); push!(pool.balance, 0); push!(pool.size, 1)
    return n
    end
end

# lc_elem_cmp on (y, i); `x` given as values.
@inline function _rmq_cmp(pool::_RMQPool, y::Int32, i::Int64, p::Int)
    py = pool.y[p]
    y < py && return -1
    y > py && return 1
    pi = pool.i[p]
    return (i > pi) - (i < pi)
end
@inline _rmq_lt2(pool::_RMQPool, a::Int, b::Int) = pool.pri[a] < pool.pri[b]
@inline _rmq_size(pool::_RMQPool, p::Int) = p == 0 ? 0 : pool.size[p]
@inline _rmq_size_child(pool::_RMQPool, q::Int, d::Int) = (c = pool.child[d + 1, q]; c == 0 ? 0 : pool.size[c])

function _rmq_update_min!(pool::_RMQPool, p::Int, q::Int, r::Int)
    @inbounds begin
    pool.s[p] = (q == 0 || _rmq_lt2(pool, p, pool.s[q])) ? p : pool.s[q]
    pool.s[p] = (r == 0 || _rmq_lt2(pool, pool.s[p], pool.s[r])) ? pool.s[p] : pool.s[r]
    return nothing
    end
end

function _rmq_rotate1!(pool::_RMQPool, p::Int, dir::Int)
    @inbounds begin
    opp = 1 - dir
    q = pool.child[opp + 1, p]; s = pool.s[p]
    size_p = pool.size[p]
    pool.size[p] -= pool.size[q] - _rmq_size_child(pool, q, dir)
    pool.size[q] = size_p
    _rmq_update_min!(pool, p, pool.child[dir + 1, p], pool.child[dir + 1, q])
    pool.s[q] = s
    pool.child[opp + 1, p] = pool.child[dir + 1, q]
    pool.child[dir + 1, q] = p
    return q
    end
end

function _rmq_rotate2!(pool::_RMQPool, p::Int, dir::Int)
    @inbounds begin
    opp = 1 - dir
    q = pool.child[opp + 1, p]; r = pool.child[dir + 1, q]; s = pool.s[p]
    size_x_dir = _rmq_size_child(pool, r, dir)
    pool.size[r] = pool.size[p]
    pool.size[p] -= pool.size[q] - size_x_dir
    pool.size[q] -= size_x_dir + 1
    _rmq_update_min!(pool, p, pool.child[dir + 1, p], pool.child[dir + 1, r])
    _rmq_update_min!(pool, q, pool.child[opp + 1, q], pool.child[opp + 1, r])
    pool.s[r] = s
    pool.child[opp + 1, p] = pool.child[dir + 1, r]
    pool.child[dir + 1, r] = p
    pool.child[dir + 1, q] = pool.child[opp + 1, r]
    pool.child[opp + 1, r] = q
    b1 = dir == 0 ? Int8(1) : Int8(-1)
    if pool.balance[r] == b1
        pool.balance[q] = 0; pool.balance[p] = -b1
    elseif pool.balance[r] == 0
        pool.balance[q] = 0; pool.balance[p] = 0
    else
        pool.balance[q] = b1; pool.balance[p] = 0
    end
    pool.balance[r] = 0
    return r
    end
end

# krmq_insert; `root` is a one-element vector holding the root index.
function _rmq_insert!(pool::_RMQPool, root::Base.RefValue{Int}, x::Int)
    @inbounds begin
    stack = pool.stack_insert; path = pool.path
    bp = root[]; bq = 0
    p = bp; q = 0; top = 0; path_len = 0; which = 0
    y = pool.y[x]; i = pool.i[x]
    while p != 0
        cmp_result = _rmq_cmp(pool, y, i, p)
        cmp_result == 0 && return p
        if pool.balance[p] != 0
            bq = q; bp = p; top = 0
        end
        which = cmp_result > 0 ? 1 : 0
        top += 1; stack[top] = which
        path_len += 1; path[path_len] = p
        q = p; p = pool.child[which + 1, p]
    end
    pool.balance[x] = 0; pool.size[x] = 1; pool.child[1, x] = 0; pool.child[2, x] = 0; pool.s[x] = x
    if q == 0
        root[] = x
    else
        pool.child[which + 1, q] = x
    end
    bp == 0 && return x
    for k in 1:path_len
        pool.size[path[k]] += 1
    end
    for k in path_len:-1:1
        _rmq_update_min!(pool, path[k], pool.child[1, path[k]], pool.child[2, path[k]])
        pool.s[path[k]] != x && break
    end
    p = bp; top = 1
    while p != x
        if stack[top] == 0
            pool.balance[p] -= 1
        else
            pool.balance[p] += 1
        end
        p = pool.child[stack[top] + 1, p]; top += 1
    end
    -2 < pool.balance[bp] < 2 && return x
    which = pool.balance[bp] < 0 ? 1 : 0
    b1 = which == 0 ? Int8(1) : Int8(-1)
    q = pool.child[2 - which, bp]
    if pool.balance[q] == b1
        r = _rmq_rotate1!(pool, bp, which)
        pool.balance[q] = 0; pool.balance[bp] = 0
    else
        r = _rmq_rotate2!(pool, bp, which)
    end
    if bq == 0
        root[] = r
    else
        pool.child[(bp != pool.child[1, bq] ? 1 : 0) + 1, bq] = r
    end
    return x
    end
end

function _rmq_find(pool::_RMQPool, root::Int, y::Int32, i::Int64)
    @inbounds begin
    p = root
    while p != 0
        c = _rmq_cmp(pool, y, i, p)
        c < 0 ? (p = pool.child[1, p]) : c > 0 ? (p = pool.child[2, p]) : break
    end
    return p
    end
end

# krmq_erase of node x (present in the tree). As in C, the walk starts from a
# fake node (pool node 1) whose left child is the root.
function _rmq_erase!(pool::_RMQPool, root::Base.RefValue{Int}, x::Int)
    @inbounds begin
    fake = 1
    pool.child[1, fake] = root[]; pool.child[2, fake] = 0
    path = pool.path; dir = pool.dir
    d = 0
    p = fake; cmp_result = -1
    y = pool.y[x]; i = pool.i[x]
    while cmp_result != 0
        which = cmp_result > 0 ? 1 : 0
        d += 1; dir[d] = which; path[d] = p
        p = pool.child[which + 1, p]
        p == 0 && return 0
        cmp_result = _rmq_cmp(pool, y, i, p)
    end
    for k in 2:d
        pool.size[path[k]] -= 1
    end
    if pool.child[2, p] == 0
        pool.child[dir[d] + 1, path[d]] = pool.child[1, p]
    else
        q = pool.child[2, p]
        if pool.child[1, q] == 0
            pool.child[1, q] = pool.child[1, p]
            pool.balance[q] = pool.balance[p]
            pool.child[dir[d] + 1, path[d]] = q
            d += 1; path[d] = q; dir[d] = 1
            pool.size[q] = pool.size[p] - 1
        else
            e = d + 1; d += 1       # slot e is reserved for the successor r
            local r
            while true
                d += 1; dir[d] = 0; path[d] = q
                r = pool.child[1, q]
                pool.child[1, r] == 0 && break
                q = r
            end
            pool.child[1, r] = pool.child[1, p]
            pool.child[1, q] = pool.child[2, r]
            pool.child[2, r] = pool.child[2, p]
            pool.balance[r] = pool.balance[p]
            pool.child[dir[e - 1] + 1, path[e - 1]] = r
            path[e] = r; dir[e] = 1
            for k in (e + 1):d
                pool.size[path[k]] -= 1
            end
            pool.size[r] = pool.size[p] - 1
        end
    end
    for k in d:-1:2       # the fake node's minimum is never read
        n = path[k]
        _rmq_update_min!(pool, n, pool.child[1, n], pool.child[2, n])
    end
    while d > 1
        q = path[d]
        which = dir[d]; other = 1 - which
        b1 = Int8(1); b2 = Int8(2)
        if which == 1
            b1 = -b1; b2 = -b2
        end
        pool.balance[q] += b1
        if pool.balance[q] == b1
            break
        elseif pool.balance[q] == b2
            r = pool.child[other + 1, q]
            if pool.balance[r] == -b1
                pool.child[dir[d - 1] + 1, path[d - 1]] = _rmq_rotate2!(pool, q, which)
            else
                pool.child[dir[d - 1] + 1, path[d - 1]] = _rmq_rotate1!(pool, q, which)
                if pool.balance[r] == 0
                    pool.balance[r] = -b1
                    pool.balance[q] = b1
                    break
                else
                    pool.balance[r] = 0; pool.balance[q] = 0
                end
            end
        end
        d -= 1
    end
    root[] = pool.child[1, fake]
    return p
    end
end

# krmq_rmq on the closed interval [lo, hi] of keys; returns 0 if empty.
function _rmq_query(pool::_RMQPool, root::Int, lo_y::Int32, lo_i::Int64, hi_y::Int32, hi_i::Int64)
    @inbounds begin
    root == 0 && return 0
    path0 = pool.path; cmp0 = pool.cmp1; path1 = pool.path2; cmp1 = pool.cmp2
    n0 = 0; n1 = 0
    p = root
    while p != 0
        c = _rmq_cmp(pool, lo_y, lo_i, p)
        n0 += 1; path0[n0] = p; cmp0[n0] = c
        c < 0 ? (p = pool.child[1, p]) : c > 0 ? (p = pool.child[2, p]) : break
    end
    p = root
    while p != 0
        c = _rmq_cmp(pool, hi_y, hi_i, p)
        n1 += 1; path1[n1] = p; cmp1[n1] = c
        c < 0 ? (p = pool.child[1, p]) : c > 0 ? (p = pool.child[2, p]) : break
    end
    k = 1
    while k <= n0 && k <= n1
        (path0[k] == path1[k] && cmp0[k] <= 0 && cmp1[k] >= 0) && break
        k += 1
    end
    (k > n0 || k > n1) && return 0
    lca = k
    best = path0[lca]
    for k in (lca + 1):n0
        if cmp0[k] <= 0
            _rmq_lt2(pool, path0[k], best) && (best = path0[k])
            right = pool.child[2, path0[k]]
            (right != 0 && _rmq_lt2(pool, pool.s[right], best)) && (best = pool.s[right])
        end
    end
    for k in (lca + 1):n1
        if cmp1[k] >= 0
            _rmq_lt2(pool, path1[k], best) && (best = path1[k])
            left = pool.child[1, path1[k]]
            (left != 0 && _rmq_lt2(pool, pool.s[left], best)) && (best = pool.s[left])
        end
    end
    return best
    end
end

# krmq_interval: (lower, upper) neighbours of the key.
function _rmq_interval(pool::_RMQPool, root::Int, y::Int32, i::Int64)
    @inbounds begin
    p = root; l = 0; u = 0
    while p != 0
        c = _rmq_cmp(pool, y, i, p)
        if c < 0
            u = p; p = pool.child[1, p]
        elseif c > 0
            l = p; p = pool.child[2, p]
        else
            l = u = p
            break
        end
    end
    return l, u
    end
end

# Iterator stack (krmq_itr_find / krmq_itr_next_bidir).
function _rmq_itr_find(pool::_RMQPool, root::Int, x::Int)
    @inbounds begin
    stack = pool.stack
    empty!(stack)
    p = root
    y = pool.y[x]; i = pool.i[x]
    while p != 0
        push!(stack, p)
        c = _rmq_cmp(pool, y, i, p)
        c < 0 ? (p = pool.child[1, p]) : c > 0 ? (p = pool.child[2, p]) : break
    end
    return stack
    end
end

function _rmq_itr_prev!(pool::_RMQPool, stack::Vector{Int})
    @inbounds begin
    isempty(stack) && return false
    p = pool.child[1, stack[end]]
    if p != 0
        while p != 0
            push!(stack, p)
            p = pool.child[2, p]
        end
        return true
    else
        while true
            p = pop!(stack)
            (!isempty(stack) && p == pool.child[1, stack[end]]) || break
        end
        return !isempty(stack)
    end
    end
end
