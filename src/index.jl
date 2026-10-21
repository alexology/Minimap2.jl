# In-memory minimizer index (index.c, `mm_idx_str` and lookups).

"""
    MinimapIndex(sequences; names=nothing, k=15, w=10, hpc=false, bucket_bits=14)
    MinimapIndex(sequences, options::IndexOptions; names=nothing)

Minimizer index of the reference `sequences` (byte vectors or strings), as built
by minimap2's `mm_idx_str`. Reference sequences are stored as 2-bit codes with 4
for ambiguous bases.
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
    table::Dict{UInt64,Vector{UInt64}}
end

function MinimapIndex(sequences::AbstractVector; names=nothing, k::Integer=15, w::Integer=10,
                      hpc::Bool=false, bucket_bits::Integer=14)
    isempty(sequences) && throw(ArgumentError("no reference sequence"))
    b = k * 2 < bucket_bits ? 2k : bucket_bits
    w = max(w, 1)
    names === nothing || length(names) == length(sequences) || throw(ArgumentError("one name per sequence"))
    entries = MM128[]
    codes = UInt8[]; lengths = Int[]; offsets = Int[]
    for (rid, sequence) in enumerate(sequences)
        bytes = sequence isa AbstractString ? codeunits(sequence) : sequence
        push!(offsets, length(codes)); push!(lengths, length(bytes))
        append!(codes, nt4.(bytes))
        isempty(bytes) || sketch!(entries, bytes, w, k, rid - 1, hpc)
    end
    table = Dict{UInt64,Vector{UInt64}}()
    for entry in entries
        push!(get!(Vector{UInt64}, table, entry.x >> 8), entry.y)
    end
    foreach(sort!, values(table))
    seq_names = names === nothing ? ["" for _ in sequences] : String.(names)
    return MinimapIndex(w, k, b, hpc, seq_names, lengths, offsets, codes, table)
end
MinimapIndex(sequences::AbstractVector, options::IndexOptions; names=nothing) =
    MinimapIndex(sequences; names, k=options.k, w=options.w, hpc=options.hpc, bucket_bits=options.bucket_bits)

const _NO_HITS = UInt64[]
# mm_idx_get: positions (rid << 32 | last_position << 1 | strand) of a minimizer hash.
@inline index_get(index::MinimapIndex, minimizer::UInt64) = get(index.table, minimizer, _NO_HITS)

"""
    max_occurrence(index, fraction) -> Int

Occurrence threshold above which the top `fraction` of distinct minimizers lie
(minimap2 `mm_idx_cal_max_occ`).
"""
function max_occurrence(index::MinimapIndex, fraction::Float32)
    fraction <= 0 && return Int(typemax(Int32))
    counts = UInt32[length(positions) for positions in values(index.table)]
    isempty(counts) && return Int(typemax(Int32))
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
