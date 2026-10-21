# Store minimap2 outputs in test/reference.
#   julia tools/make_reference.jl [minimap2 executable]
# Default executable: reference/minimap2/minimap2(.exe), built from the pinned
# checkout (see README).
using Pkg
Pkg.activate(; temp=true, io=devnull); Pkg.add("CodecZlib"; io=devnull)
using CodecZlib

include(joinpath(@__DIR__, "..", "test", "cases.jl"))

function main()
    minimap2 = isempty(ARGS) ? joinpath(@__DIR__, "..", "reference", "minimap2", Sys.iswindows() ? "minimap2.exe" : "minimap2") : ARGS[1]
    mkpath(REFERENCE)
    store(name, bytes) = (write(joinpath(REFERENCE, name), transcode(GzipCompressor, bytes)); println(name))
    for (dataset, refs, reads) in DATASETS, preset in PRESETS
        for cigar in (false, true)
            out = read(pipeline(paf_command(minimap2, refs, reads, preset, cigar); stderr=devnull), String)
            store(paf_name(dataset, preset, cigar), Vector{UInt8}(replace(out, "\r" => "")))
        end
        preset == "map-ont" || continue    # seeds and chains for one preset only (size)
        err = IOBuffer()
        run(pipeline(debug_command(minimap2, refs, reads, preset); stdout=devnull, stderr=err))
        lines = filter(l -> startswith(l, "QR\t") || startswith(l, "RS\t") || startswith(l, "SD\t") || startswith(l, "CN\t"),
                       split(replace(String(take!(err)), "\r" => ""), '\n'))
        store(debug_name(dataset, preset), Vector{UInt8}(join(lines, '\n') * "\n"))
    end
end
main()
