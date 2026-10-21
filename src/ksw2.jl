# ksw2 alignment kernels used by minimap2, emulated cell by cell with the
# arithmetic of their SSE implementations (8-bit wrapping differences in
# ksw_extd2_sse, 16-bit saturating scores in ksw_ll_i16), so that scores,
# coordinates and CIGARs, including tie-breaks, are identical.

const KSW_NEG_INF = -Int32(0x40000000)
const KSW_EZ_SCORE_ONLY = Int(0x01)
const KSW_EZ_RIGHT = Int(0x02)
const KSW_EZ_GENERIC_SC = Int(0x04)
const KSW_EZ_APPROX_MAX = Int(0x08)
const KSW_EZ_APPROX_DROP = Int(0x10)
const KSW_EZ_EXTZ_ONLY = Int(0x40)
const KSW_EZ_REV_CIGAR = Int(0x80)

const CIGAR_MATCH = 0
const CIGAR_INS = 1
const CIGAR_DEL = 2
const CIGAR_N_SKIP = 3

"""
    KswResult

Result of an extension alignment (ksw2 `ksw_extz_t`).
"""
mutable struct KswResult
    max::Int32
    zdropped::Bool
    max_q::Int32
    max_t::Int32
    mqe::Int32
    mqe_t::Int32
    mte::Int32
    mte_q::Int32
    score::Int32
    reach_end::Bool
    cigar::Vector{UInt32}
    # Work buffers reused across calls.
    work8::Vector{Vector{Int8}}
    sf::Vector{UInt8}
    qr::Vector{UInt8}
    H::Vector{Int32}
    p::Vector{UInt8}
    off::Vector{Int}
    off_end::Vector{Int}
end
KswResult() = KswResult(0, false, -1, -1, KSW_NEG_INF, -1, KSW_NEG_INF, -1, KSW_NEG_INF, false, UInt32[],
                        [Int8[] for _ in 1:10], UInt8[], UInt8[], Int32[], UInt8[], Int[], Int[])

function _reset!(ez::KswResult)
    ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1
    ez.max = 0; ez.score = ez.mqe = ez.mte = KSW_NEG_INF
    empty!(ez.cigar); ez.zdropped = false; ez.reach_end = false
    return ez
end

function _push_cigar!(cigar::Vector{UInt32}, op::Integer, len::Integer)
    if isempty(cigar) || op != (cigar[end] & 0xf)
        push!(cigar, UInt32(len) << 4 | UInt32(op))
    else
        cigar[end] += UInt32(len) << 4
    end
    return cigar
end

# ksw_apply_zdrop with is_rot = 1 (a = anti-diagonal, b = target position).
function _apply_zdrop!(ez::KswResult, H::Int32, a::Int, b::Int, zdrop::Int, e::Int)
    r = a; t = b
    if H > ez.max
        ez.max = H; ez.max_t = t; ez.max_q = r - t
    elseif t >= ez.max_t && r - t >= ez.max_q
        tl = t - ez.max_t; ql = (r - t) - ez.max_q
        l = tl > ql ? tl - ql : ql - tl
        if zdrop >= 0 && ez.max - H > zdrop + l * e
            ez.zdropped = true
            return true
        end
    end
    return false
end

# ksw_backtrack for the rotated (anti-diagonal) layout; i = target, j = query.
function _ksw_backtrack(is_rev::Bool, p::Vector{UInt8}, n_col::Int, off::Vector{Int}, off_end::Vector{Int}, i0::Int, j0::Int)
    cigar = UInt32[]
    i = i0; j = j0; state = 0
    while i >= 0 && j >= 0
        force_state = -1
        r = i + j
        i < off[r + 1] && (force_state = 2)
        i > off_end[r + 1] && (force_state = 1)
        tmp = force_state < 0 ? Int(p[r * n_col + i - off[r + 1] + 1]) : 0
        if state == 0
            state = tmp & 7
        elseif ((tmp >> (state + 2)) & 1) == 0
            state = 0
        end
        state == 0 && (state = tmp & 7)
        force_state >= 0 && (state = force_state)
        if state == 0
            _push_cigar!(cigar, CIGAR_MATCH, 1); i -= 1; j -= 1
        elseif state == 1 || state == 3
            _push_cigar!(cigar, CIGAR_DEL, 1); i -= 1
        else
            _push_cigar!(cigar, CIGAR_INS, 1); j -= 1
        end
    end
    i >= 0 && _push_cigar!(cigar, CIGAR_DEL, i + 1)
    j >= 0 && _push_cigar!(cigar, CIGAR_INS, j + 1)
    is_rev || reverse!(cigar)
    return cigar
end
