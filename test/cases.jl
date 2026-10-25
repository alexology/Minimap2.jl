# Datasets and presets compared with minimap2 2.30. Shared by test/runtests.jl and
# tools/make_reference.jl.

const DATA = joinpath(@__DIR__, "data")
const REFERENCE = joinpath(@__DIR__, "reference")
# (name, references, reads)
const DATASETS = [("amplicon", "amplicon_refs.fa", "amplicon_reads.fq"),
                  ("noisy", "noisy_refs.fa", "noisy_reads.fq"),
                  ("hard", "hard_refs.fa", "hard_reads.fa"),
                  # the inverted part ends in the padding of the striped query profile
                  ("inversion", "inversion_refs.fa", "inversion_reads.fa")]
const PRESETS = ["map-ont", "lr:hq", "map-pb"]

paf_name(dataset, preset, cigar) = "$(dataset)_$(replace(preset, ':' => '_'))$(cigar ? "_cigar" : "").paf.gz"
debug_name(dataset, preset) = "$(dataset)_$(replace(preset, ':' => '_'))_debug.tsv.gz"
paf_command(minimap2, refs, reads, preset, cigar) =
    `$minimap2 -x $preset -t 1 $(cigar ? ["-c"] : String[]) $(joinpath(DATA, refs)) $(joinpath(DATA, reads))`
debug_command(minimap2, refs, reads, preset) =
    `$minimap2 -x $preset -t 1 --print-seeds --print-chains $(joinpath(DATA, refs)) $(joinpath(DATA, reads))`
