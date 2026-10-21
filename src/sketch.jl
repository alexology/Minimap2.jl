# Symmetric (w,k)-minimizers (sketch.c).

const SEQ_NT4 = let table = fill(UInt8(4), 256)
    for (base, code) in zip("ACGTUacgtu", (0, 1, 2, 3, 3, 0, 1, 2, 3, 3))
        table[UInt8(base) + 1] = code
    end
    table
end
@inline nt4(base::UInt8) = @inbounds SEQ_NT4[base + 1]

@inline function _hash64(key::UInt64, mask::UInt64)
    key = (~key + (key << 21)) & mask
    key = key ⊻ (key >> 24)
    key = ((key + (key << 3)) + (key << 8)) & mask
    key = key ⊻ (key >> 14)
    key = ((key + (key << 2)) + (key << 4)) & mask
    key = key ⊻ (key >> 28)
    key = (key + (key << 31)) & mask
    return key
end

"""
    sketch!(minimizers, sequence, w, k, rid, hpc) -> minimizers

Append the symmetric (w,k)-minimizers of `sequence` to `minimizers` (minimap2
`mm_sketch`). Each entry has `x = hash << 8 | span` and
`y = rid << 32 | last_position << 1 | strand`, with 0-based positions. With
`hpc`, k-mers are taken on the homopolymer-compressed sequence.
"""
function sketch!(minimizers::Vector{MM128}, sequence::AbstractVector{UInt8}, w::Integer, k::Integer,
                 rid::Integer, hpc::Bool)
    len = length(sequence)
    (len > 0 && 0 < w < 256 && 0 < k <= 28) || throw(ArgumentError("invalid sketch parameters"))
    shift1 = 2 * (k - 1); mask = (UInt64(1) << (2k)) - 1
    forward = UInt64(0); reverse = UInt64(0)
    none = MM128(typemax(UInt64), typemax(UInt64))
    buffer = fill(none, w)
    minimum_entry = none
    queue = zeros(Int, 32); queue_front = 0; queue_count = 0
    l = 0; buf_pos = 0; min_pos = 0; kmer_span = 0
    i = 0
    while i < len
        c = nt4(sequence[i + 1])
        info = none
        if c < 4
            if hpc
                skip = 1
                if i + 1 < len && nt4(sequence[i + 2]) == c
                    skip = 2
                    while i + skip < len
                        nt4(sequence[i + skip + 1]) != c && break
                        skip += 1
                    end
                    i += skip - 1
                end
                queue[((queue_count + queue_front) & 0x1f) + 1] = skip
                queue_count += 1
                kmer_span += skip
                if queue_count > k
                    kmer_span -= queue[queue_front + 1]
                    queue_front = (queue_front + 1) & 0x1f
                    queue_count -= 1
                end
            else
                kmer_span = l + 1 < k ? l + 1 : k
            end
            forward = ((forward << 2) | c) & mask
            reverse = (reverse >> 2) | (UInt64(3 ⊻ c) << shift1)
            if forward == reverse
                i += 1
                continue
            end
            z = forward < reverse ? 0 : 1
            l += 1
            if l >= k && kmer_span < 256
                info = MM128(_hash64(z == 0 ? forward : reverse, mask) << 8 | UInt64(kmer_span),
                             UInt64(rid) << 32 | UInt64(UInt32(i)) << 1 | UInt64(z))
            end
        else
            l = 0; queue_count = 0; queue_front = 0; kmer_span = 0
        end
        buffer[buf_pos + 1] = info
        if l == w + k - 1 && minimum_entry.x != typemax(UInt64)
            for j in (buf_pos + 1):(w - 1)
                (minimum_entry.x == buffer[j + 1].x && buffer[j + 1].y != minimum_entry.y) && push!(minimizers, buffer[j + 1])
            end
            for j in 0:(buf_pos - 1)
                (minimum_entry.x == buffer[j + 1].x && buffer[j + 1].y != minimum_entry.y) && push!(minimizers, buffer[j + 1])
            end
        end
        if info.x <= minimum_entry.x
            (l >= w + k && minimum_entry.x != typemax(UInt64)) && push!(minimizers, minimum_entry)
            minimum_entry = info; min_pos = buf_pos
        elseif buf_pos == min_pos
            (l >= w + k - 1 && minimum_entry.x != typemax(UInt64)) && push!(minimizers, minimum_entry)
            minimum_entry = none
            for j in (buf_pos + 1):(w - 1)
                if minimum_entry.x >= buffer[j + 1].x
                    minimum_entry = buffer[j + 1]; min_pos = j
                end
            end
            for j in 0:buf_pos
                if minimum_entry.x >= buffer[j + 1].x
                    minimum_entry = buffer[j + 1]; min_pos = j
                end
            end
            if l >= w + k - 1 && minimum_entry.x != typemax(UInt64)
                for j in (buf_pos + 1):(w - 1)
                    (minimum_entry.x == buffer[j + 1].x && minimum_entry.y != buffer[j + 1].y) && push!(minimizers, buffer[j + 1])
                end
                for j in 0:buf_pos
                    (minimum_entry.x == buffer[j + 1].x && minimum_entry.y != buffer[j + 1].y) && push!(minimizers, buffer[j + 1])
                end
            end
        end
        buf_pos += 1
        buf_pos == w && (buf_pos = 0)
        i += 1
    end
    minimum_entry.x != typemax(UInt64) && push!(minimizers, minimum_entry)
    return minimizers
end
