# In-memory minimizer index (index.c, `mm_idx_str` and lookups).

"""
    MinimapIndex(sequences; names=nothing, k=15, w=10, hpc=false, bucket_bits=14, threads=Threads.nthreads())
    MinimapIndex(sequences, options::IndexOptions; names=nothing, threads=Threads.nthreads())

Minimizer index of the reference `sequences` (byte vectors or strings), as built
by minimap2's `mm_idx_str`. Reference sequences are stored as 2-bit codes with 4
for ambiguous bases.

Minimizers are stored sorted by hash in flat arrays: the positions of the `i`-th
distinct hash `keys[i]` are `positions[starts[i]+1:starts[i+1]]`, in increasing
order, and `directory` maps the high bits of a hash to its range of `keys`.
`threads` tasks sketch and sort the references; the index does not depend on it.
"""
struct MinimapIndex
    w::Int
    k::Int
    b::Int
    hpc::Bool
    names::Vector{String}
    lengths::Vector{Int}
    offsets::Vector{Int}
    codes::Vector{UInt8}
    keys::Vector{UInt64}
    starts::Vector{UInt32}
    positions::Vector{UInt64}
    directory::Vector{UInt32}
    directory_shift::Int
end

const _PositionView = SubArray{UInt64,1,Vector{UInt64},Tuple{UnitRange{Int}},true}

# Below this many reference bases the index is built on the calling task.
const _PARALLEL_MIN_BASES = 1 << 20
# Entries are scattered into 2^_SORT_BITS buckets of hash high bits, sorted one by one.
const _SORT_BITS = 12

function MinimapIndex(sequences::AbstractVector; names=nothing, k::Integer=15, w::Integer=10,
                      hpc::Bool=false, bucket_bits::Integer=14, threads::Integer=Threads.nthreads())
    isempty(sequences) && throw(ArgumentError("no reference sequence"))
    threads >= 1 || throw(ArgumentError("threads must be positive"))
    (0 < k <= 28) || throw(ArgumentError("invalid sketch parameters"))
    b = k * 2 < bucket_bits ? 2k : bucket_bits
    w = max(w, 1)
    names === nothing || length(names) == length(sequences) || throw(ArgumentError("one name per sequence"))
    bytes = [sequence isa AbstractString ? codeunits(sequence) : sequence for sequence in sequences]
    lengths = Int[length(s) for s in bytes]
    offsets = cumsum([0; lengths[1:(end - 1)]])
    codes = Vector{UInt8}(undef, sum(lengths))
    chunks = _reference_chunks(lengths, sum(lengths) < _PARALLEL_MIN_BASES ? 1 : threads)
    # Sketch each chunk of references and store their 2-bit codes.
    chunk_entries = _tmap(threads, chunks) do range
        entries = MM128[]
        for rid in range
            s = bytes[rid]
            for j in eachindex(s)
                @inbounds codes[offsets[rid] + j] = nt4(s[j])
            end
            isempty(s) || sketch!(entries, s, w, k, rid - 1, hpc)
        end
        entries
    end
    keys, starts, positions = _sort_entries(chunk_entries, 2k, threads)
    directory, directory_shift = _key_directory(keys, 2k)
    seq_names = names === nothing ? ["" for _ in sequences] : String.(names)
    return MinimapIndex(w, k, b, hpc, seq_names, lengths, offsets, codes, keys, starts, positions,
                        directory, directory_shift)
end
MinimapIndex(sequences::AbstractVector, options::IndexOptions; names=nothing, threads::Integer=Threads.nthreads()) =
    MinimapIndex(sequences; names, k=options.k, w=options.w, hpc=options.hpc, bucket_bits=options.bucket_bits, threads)

# Results of `f` over `items`, computed by up to `threads` tasks (in order).
function _tmap(f, threads::Integer, items::AbstractVector)
    n = length(items)
    (threads == 1 || n <= 1) && return map(f, items)
    results = Vector{Any}(undef, n)
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:min(threads, n)
        Threads.@spawn while true
            i = Threads.atomic_add!(next, 1)
            i > n && break
            results[i] = f(items[i])
        end
    end
    return [r for r in results]
end

# Contiguous reference ranges of similar total length, a few per thread.
function _reference_chunks(lengths::Vector{Int}, threads::Integer)
    n = length(lengths)
    threads == 1 && return [1:n]
    target = max(sum(lengths) ÷ (4threads), 1)
    chunks = UnitRange{Int}[]; first = 1; acc = 0
    for i in 1:n
        acc += lengths[i]
        if acc >= target || i == n
            push!(chunks, first:i); first = i + 1; acc = 0
        end
    end
    return chunks
end

# Flatten the sketches into (key, position) pairs sorted by key and position:
# distinct keys, the offset of each key's positions, and the positions.
function _sort_entries(chunk_entries::Vector{Vector{MM128}}, key_bits::Int, threads::Integer)
    sort_bits = min(_SORT_BITS, key_bits)
    shift = key_bits - sort_bits; nb = 1 << sort_bits
    nc = length(chunk_entries)
    # Counting sort on the high bits: chunk c writes bucket j at slots[j, c].
    counts = zeros(Int, nb, nc)
    _tforeach(threads, 1:nc) do c
        for e in chunk_entries[c]
            @inbounds counts[(e.x >> 8) >> shift + 1, c] += 1
        end
    end
    slots = similar(counts); bucket_start = zeros(Int, nb + 1); acc = 0
    for j in 1:nb
        bucket_start[j] = acc
        for c in 1:nc
            slots[j, c] = acc; acc += counts[j, c]
        end
    end
    bucket_start[nb + 1] = acc
    pairs = Vector{UInt128}(undef, acc)
    _tforeach(threads, 1:nc) do c
        for e in chunk_entries[c]
            j = (e.x >> 8) >> shift + 1
            @inbounds pairs[slots[j, c] += 1] = UInt128(e.x >> 8) << 64 | e.y
        end
        empty!(chunk_entries[c]); sizehint!(chunk_entries[c], 0)
    end
    # Sort each bucket and count its distinct keys.
    bucket_keys = zeros(Int, nb)
    _tforeach(threads, 1:nb) do j
        range = (bucket_start[j] + 1):bucket_start[j + 1]
        isempty(range) && return
        sort!(view(pairs, range))
        n = 1
        for i in (first(range) + 1):last(range)
            @inbounds (pairs[i] >> 64) != (pairs[i - 1] >> 64) && (n += 1)
        end
        bucket_keys[j] = n
    end
    acc <= typemax(UInt32) || throw(ArgumentError("more than $(typemax(UInt32)) minimizers"))
    key_start = cumsum([0; bucket_keys])
    keys = Vector{UInt64}(undef, key_start[end])
    starts = Vector{UInt32}(undef, key_start[end] + 1)
    positions = Vector{UInt64}(undef, acc)
    _tforeach(threads, 1:nb) do j
        range = (bucket_start[j] + 1):bucket_start[j + 1]
        isempty(range) && return
        ki = key_start[j]
        for i in range
            @inbounds begin
                key = UInt64(pairs[i] >> 64)
                if i == first(range) || key != keys[ki]
                    ki += 1; keys[ki] = key; starts[ki] = i - 1
                end
                positions[i] = pairs[i] % UInt64
            end
        end
    end
    starts[end] = acc
    return keys, starts, positions
end

function _tforeach(f, threads::Integer, items::AbstractVector)
    if threads == 1 || length(items) <= 1
        foreach(f, items)
    else
        _tmap(x -> (f(x); nothing), threads, items)
    end
    return nothing
end

# Directory over the high bits of the keys, about one key per entry:
# keys[directory[h+1]+1:directory[h+2]] have `key >> shift == h`.
function _key_directory(keys::Vector{UInt64}, key_bits::Int)
    bits = clamp(floor(Int, log2(max(length(keys), 1))), 1, key_bits)
    shift = key_bits - bits
    directory = zeros(UInt32, (1 << bits) + 1)
    for key in keys
        @inbounds directory[key >> shift + 2] += 1
    end
    for h in 2:length(directory)
        @inbounds directory[h] += directory[h - 1]
    end
    return directory, shift
end

# mm_idx_get: positions (rid << 32 | last_position << 1 | strand) of a minimizer hash.
@inline function index_get(index::MinimapIndex, minimizer::UInt64)
    h = minimizer >> index.directory_shift
    h + 2 > length(index.directory) && return view(index.positions, 1:0)
    @inbounds begin
        lo = Int(index.directory[h + 1]) + 1; hi = Int(index.directory[h + 2])
        keys = index.keys
        while lo <= hi
            mid = (lo + hi) >>> 1
            key = keys[mid]
            if key < minimizer
                lo = mid + 1
            elseif key > minimizer
                hi = mid - 1
            else
                return view(index.positions, (Int(index.starts[mid]) + 1):Int(index.starts[mid + 1]))
            end
        end
    end
    return view(index.positions, 1:0)
end

"""
    max_occurrence(index, fraction) -> Int

Occurrence threshold above which the top `fraction` of distinct minimizers lie
(minimap2 `mm_idx_cal_max_occ`).
"""
function max_occurrence(index::MinimapIndex, fraction::Float32)
    fraction <= 0 && return Int(typemax(Int32))
    isempty(index.keys) && return Int(typemax(Int32))
    counts = diff(index.starts)
    sort!(counts)
    kk = unsafe_trunc(UInt32, (1.0 - Float64(fraction)) * length(counts))
    return Int(counts[kk + 1]) + 1
end

# mm_idx_getseq: 2-bit codes of reference `rid` (0-based) in [start, stop).
function reference_codes(index::MinimapIndex, rid::Integer, start::Integer, stop::Integer)
    offset = index.offsets[rid + 1]
    stop = min(stop, index.lengths[rid + 1])
    start >= stop && return UInt8[]
    return index.codes[(offset + start + 1):(offset + stop)]
end
