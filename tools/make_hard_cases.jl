# Synthetic references and reads with inversions, chimeras, long indels, N runs and
# high error rates (development check data).
#   julia tools/make_hard_cases.jl <output directory>
mutable struct Rng; s::UInt64; end
function nxt!(r::Rng)
    r.s += 0x9e3779b97f4a7c15; z = r.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9; z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end
unif(r) = (nxt!(r) >> 11) * 2.0^-53
rint(r, n) = Int(nxt!(r) % UInt64(n)) + 1
const B = b"ACGT"
comp(c) = c == UInt8('A') ? UInt8('T') : c == UInt8('T') ? UInt8('A') : c == UInt8('C') ? UInt8('G') : c == UInt8('G') ? UInt8('C') : c
revcomp(s) = reverse!(map(comp, s))
function noisy(r, s, err)
    out = UInt8[]
    for c in s
        x = unif(r)
        if x < err / 3
            continue
        elseif x < 2err / 3
            push!(out, B[rint(r, 4)], c)
        elseif x < err
            push!(out, B[rint(r, 4)])
        else
            push!(out, c)
        end
    end
    return out
end
r = Rng(0x4a7d)
refs = [[B[rint(r, 4)] for _ in 1:n] for n in (5000, 4000, 6000)]
# A repeat shared by references 1 and 3.
rep = [B[rint(r, 4)] for _ in 1:400]
refs[1][2001:2400] = rep; refs[3][1001:1400] = rep
reads = Pair{String,Vector{UInt8}}[]
for i in 1:12
    ref = refs[rint(r, 3)]; st = rint(r, length(ref) - 3000); s = ref[st:(st + 2999)]
    kind = i % 6
    if kind == 0       # inversion of 600 bp in the middle
        s = vcat(s[1:1200], revcomp(s[1201:1800]), s[1801:end])
    elseif kind == 1   # chimera with another reference
        o = refs[rint(r, 3)]; os = rint(r, length(o) - 1600)
        s = vcat(s[1:1500], o[os:(os + 1499)])
    elseif kind == 2   # 500-bp deletion and 300-bp insertion
        s = vcat(s[1:1000], s[1501:2200], [B[rint(r, 4)] for _ in 1:300], s[2201:end])
    elseif kind == 3   # run of N
        s[1400:1460] .= UInt8('N')
    elseif kind == 4   # low quality
        s = noisy(r, s, 0.15)
    end
    s = noisy(r, s, 0.03)
    unif(r) < 0.5 && (s = revcomp(s))
    push!(reads, "h$i" => s)
end
open(joinpath(ARGS[1], "hard_refs.fa"), "w") do io
    for (i, s) in enumerate(refs); write(io, ">ref$i\n", s, "\n"); end
end
open(joinpath(ARGS[1], "hard_reads.fa"), "w") do io
    for (name, s) in reads; write(io, ">", name, "\n", s, "\n"); end
end
