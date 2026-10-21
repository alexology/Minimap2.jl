# Minimal FASTA/FASTQ reader for the development checks: (name, sequence) pairs.
function read_fastx(path)
    records = Tuple{String,Vector{UInt8}}[]
    lines = readlines(path)
    i = 1
    while i <= length(lines)
        if startswith(lines[i], '>')
            name = split(lines[i][2:end])[1]; seq = UInt8[]; i += 1
            while i <= length(lines) && !startswith(lines[i], '>')
                append!(seq, codeunits(lines[i])); i += 1
            end
            push!(records, (name, seq))
        else
            push!(records, (split(lines[i][2:end])[1], Vector{UInt8}(codeunits(lines[i + 1])))); i += 4
        end
    end
    return records
end
