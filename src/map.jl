# Mapping driver (map.c, `mm_map_frag_core`), single-segment queries.

function _query_minimizers(index::MinimapIndex, sequence::AbstractVector{UInt8}, opt::MapOptions)
    minimizers = MM128[]
    isempty(sequence) || sketch!(minimizers, sequence, index.w, index.k, 0, index.hpc)
    opt.sdust_thres > 0 && throw(ArgumentError("SDUST masking is not implemented"))
    opt.q_occ_frac > 0 && _seed_mz_filter!(minimizers, opt.mid_occ, opt.q_occ_frac)
    return minimizers
end

"""
    map_anchors(index, sequence, opt; qname=nothing) -> (anchors, rep_len)

Query minimizers, their seed hits in `index` and the anchors sorted by reference
position (the stage that `minimap2 --print-seeds` reports). `opt` must have been
updated for `index` (its `mid_occ` set).
"""
function map_anchors(index::MinimapIndex, sequence::AbstractVector{UInt8}, opt::MapOptions; qname=nothing)
    minimizers = _query_minimizers(index, sequence, opt)
    anchors, rep_len, _ = collect_seed_hits(index, opt, opt.mid_occ, minimizers, length(sequence), qname)
    return anchors, rep_len
end

# Per-read hash used to break ties between chains of equal score.
function _read_hash(qname, qlen::Int, opt::MapOptions)
    h = (qname !== nothing && (opt.flag & F_NO_HASH_NAME) == 0) ? _x31_hash(qname) : UInt32(0)
    # minimap2 converts its int seed to an unsigned 32-bit key: negative seeds wrap.
    h ⊻= _wang_hash(UInt32(qlen)) + _wang_hash(opt.seed % UInt32)
    return _wang_hash(h)
end

"""
    map_chains(index, sequence, opt; qname=nothing) -> (regions, anchors, rep_len, mini_pos)

Seeds, chaining (with minimap2's re-chaining rules) and the regions generated from
the chains, before primary/secondary selection: the stage that
`minimap2 --print-chains` reports.
"""
function map_chains(index::MinimapIndex, sequence::AbstractVector{UInt8}, opt::MapOptions; qname=nothing)
    qlen = length(sequence)
    qlen == 0 && return Region[], MM128[], 0, UInt64[]
    (opt.max_qlen > 0 && qlen > opt.max_qlen) && return Region[], MM128[], 0, UInt64[]
    is_splice = (opt.flag & F_SPLICE) != 0; is_sr = (opt.flag & F_SR) != 0
    hash = _read_hash(qname, qlen, opt)
    minimizers = _query_minimizers(index, sequence, opt)
    a, rep_len, mini_pos = collect_seed_hits(index, opt, opt.mid_occ, minimizers, qlen, qname)
    max_chain_gap_qry = is_sr ? max(qlen, opt.max_gap) : opt.max_gap
    if opt.max_gap_ref > 0
        max_chain_gap_ref = opt.max_gap_ref
    elseif opt.max_frag_len > 0
        max_chain_gap_ref = max(opt.max_frag_len - qlen, opt.max_gap)
    else
        max_chain_gap_ref = opt.max_gap
    end
    pen_gap = Float32(Float64(opt.chain_gap_scale) * 0.01 * index.k)
    pen_skip = Float32(Float64(opt.chain_skip_scale) * 0.01 * index.k)
    if (opt.flag & F_RMQ) != 0
        u, a = chain_rmq(a, opt.max_gap, opt.rmq_inner_dist, opt.bw, opt.max_chain_skip, opt.rmq_size_cap,
                         opt.min_cnt, opt.min_chain_score, pen_gap, pen_skip)
    else
        u, a = chain_dp(a, max_chain_gap_ref, max_chain_gap_qry, opt.bw, opt.max_chain_skip, opt.max_chain_iter,
                        opt.min_cnt, opt.min_chain_score, pen_gap, pen_skip, is_splice, 1)
    end
    if opt.bw_long > opt.bw && (opt.flag & (F_SPLICE | F_SR | F_NO_LJOIN)) == 0 && length(u) > 1
        # Re-chain with RMQ and the long band when the best chain covers little of a long read.
        st = Int(_ay(a[1])); en = Int(_ay(a[Int(u[1] & 0xffffffff)]))
        if qlen - (en - st) > opt.rmq_rescue_size || en - st > qlen * opt.rmq_rescue_ratio
            n_a = sum(Int(x & 0xffffffff) for x in u)
            a = a[1:n_a]
            radix_sort_128x!(a)
            u, a = chain_rmq(a, opt.max_gap, opt.rmq_inner_dist, opt.bw_long, opt.max_chain_skip, opt.rmq_size_cap,
                             opt.min_cnt, opt.min_chain_score, pen_gap, pen_skip)
        end
    elseif opt.max_occ > opt.mid_occ && rep_len > 0 && (opt.flag & F_RMQ) == 0
        rechain = isempty(u)    # with one segment the best chain always holds all segments
        if rechain
            a, rep_len, mini_pos = collect_seed_hits(index, opt, opt.max_occ, minimizers, qlen, qname)
            u, a = chain_dp(a, max_chain_gap_ref, max_chain_gap_qry, opt.bw, opt.max_chain_skip, opt.max_chain_iter,
                            opt.min_cnt, opt.min_chain_score, pen_gap, pen_skip, is_splice, 1)
        end
    end
    regions = generate_regions(hash, qlen, u, a, (opt.flag & F_QSTRAND) != 0)
    return regions, a, rep_len, mini_pos
end

"""
    map_query(index, sequence, opt; qname=nothing) -> (regions, rep_len)

Map one query (minimap2 `mm_map`): chaining, primary/secondary selection,
divergence estimate, base-level alignment with `F_CIGAR` and mapping quality.
`qname` seeds the hash that breaks ties between equal chains, as on minimap2's
command line; without it the seed is that of `--no-hash-name`.
"""
function map_query(index::MinimapIndex, sequence::AbstractVector{UInt8}, opt::MapOptions; qname=nothing)
    regions, a, rep_len, mini_pos = map_chains(index, sequence, opt; qname)
    qlen = length(sequence)
    is_sr = (opt.flag & F_SR) != 0; is_splice = (opt.flag & F_SPLICE) != 0
    is_sr_rna = (opt.flag & F_SR_RNA) != 0
    if (opt.flag & F_ALL_CHAINS) == 0
        set_parent!(regions, opt.mask_level, opt.mask_len, opt.a * 2 + opt.b, (opt.flag & F_HARD_MLEVEL) != 0, opt.alt_drop)
        select_sub!(regions, opt.pri_ratio, index.k * 2, opt.best_n, true, unsafe_trunc(Int, opt.max_gap * 0.8))
    end
    if !is_sr && (opt.flag & F_QSTRAND) == 0
        estimate_divergence!(regions, index, qlen, a, mini_pos)
        filter_strand_retained!(regions)
    end
    if (opt.flag & F_CIGAR) != 0
        regions, a = align_regions!(regions, index, sequence, a, opt)
        if (opt.flag & F_ALL_CHAINS) == 0
            set_parent!(regions, opt.mask_level, opt.mask_len, opt.a * 2 + opt.b, (opt.flag & F_HARD_MLEVEL) != 0, opt.alt_drop)
            select_sub!(regions, opt.pri_ratio, index.k * 2, opt.best_n, false, unsafe_trunc(Int, opt.max_gap * 0.8))
            set_sam_primary!(regions)
        end
    end
    set_mapq!(regions, opt.min_chain_score, opt.a, rep_len, is_sr || is_sr_rna, is_splice)
    return regions, rep_len
end

"""
    paf_line(index, qname, qlen, region, rep_len; flag=0) -> String

The PAF record minimap2 writes for `region` (`mm_write_paf3`), without trailing
newline.
"""
function paf_line(index::MinimapIndex, qname::AbstractString, qlen::Integer, r::Region, rep_len::Integer; flag::Int64=Int64(0))
    io = IOBuffer()
    name = index.names[r.rid + 1]
    print(io, qname, '\t', qlen, '\t', r.qs, '\t', r.qe, '\t', r.rev ? '-' : '+', '\t',
          isempty(name) ? string(r.rid) : name, '\t', index.lengths[r.rid + 1])
    if (flag & F_QSTRAND) != 0 && r.rev
        print(io, '\t', index.lengths[r.rid + 1] - r.re, '\t', index.lengths[r.rid + 1] - r.rs)
    else
        print(io, '\t', r.rs, '\t', r.re)
    end
    print(io, '\t', r.mlen, '\t', r.blen, '\t', r.mapq)
    type = r.id == r.parent ? (r.inv ? 'I' : 'P') : (r.inv ? 'i' : 'S')
    p = r.extra
    if p !== nothing
        print(io, "\tNM:i:", r.blen - r.mlen + p.n_ambi, "\tms:i:", p.dp_max0, "\tAS:i:", p.dp_score, "\tnn:i:", p.n_ambi)
        (p.trans_strand == 1 || p.trans_strand == 2) && print(io, "\tts:A:", "?+-?"[p.trans_strand + 1])
    end
    print(io, "\ttp:A:", type, "\tcm:i:", r.cnt, "\ts1:i:", r.score)
    r.parent == r.id && print(io, "\ts2:i:", r.subsc)
    if p !== nothing
        div = 1.0 - event_identity(r)
        print(io, "\tde:f:", div == 0 ? "0" : @sprintf("%.4f", div))
    elseif r.div >= 0 && r.div <= 1
        print(io, "\tdv:f:", r.div == 0 ? "0" : @sprintf("%.4f", Float64(r.div)))
    end
    r.split != 0 && print(io, "\tzd:i:", r.split)
    rep_len >= 0 && print(io, "\trl:i:", rep_len)
    if p !== nothing && (flag & F_OUT_CG) != 0
        print(io, "\tcg:Z:")
        for c in p.cigar
            print(io, c >> 4, "MIDNSHP=XB"[Int(c & 0xf) + 1])
        end
    end
    return String(take!(io))
end
