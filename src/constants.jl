# constants

# log2(e)
const LOG2E = 1.442695040888963407359924681001892137426646
# log2(10)
const LOG210 = 3.321928094887362347870319429489390175864831393024580612054756395815934776608624
# log(10)
const LN10 = 2.302585092994045684017991454684364207601101488628772976033327900967572609677367


# log(2)
LN2U{T}(::Type{T}) = T(6.93145751953125e-1)
LN2L{T}(::Type{T}) = T(1.42860682030941723212e-6)

LN2U{T<:SmallFloatTypes}(::Type{T}) = T(0.693359375)
LN2L{T<:SmallFloatTypes}(::Type{T}) = T(-2.12194440e-4)

# log10(2)
LOG102U{T}(::Type{T}) = T(3.01025390625000000000e-1)
LOG102L{T}(::Type{T}) = T(4.60503898119521373889e-6)

LOG102U{T<:SmallFloatTypes}(::Type{T}) = T(3.00781250000000000000e-1)
LOG102L{T<:SmallFloatTypes}(::Type{T}) = T(2.48745663981195213739e-4)

# max and min arguments for exponential fucntions
MAXEXP(::Type{Float64}) = 7.09782712893383996732e2 # log 2^1023*(2-2^-52)
MAXEXP(::Type{Float32}) = 88.72283905206835f0      # log 2^127 *(2-2^-23)
MAXEXP(::Type{Float16}) = Float16(11.09)           # log 2^15  *(2-2^-10)

MAXEXP2(::Type{Float64}) = 1024         # log2 2^1023*(2-2^-52)
MAXEXP2(::Type{Float32}) = 128f0        # log2 2^127 *(2-2^-23)
MAXEXP2(::Type{Float16}) = Float16(16)  # log2 2^15  *(2-2^-10)

MAXEXP10(::Type{Float64}) = 3.08254715559916743851e2 # log 2^1023*(2-2^-52)
MAXEXP10(::Type{Float32}) = 38.531839419103626f0     # log 2^127 *(2-2^-23)
MAXEXP10(::Type{Float16}) = Float16(4.816268)        # log 2^15  *(2-2^-10)

# one less than the min exponent since we can sqeeze a bit more from the exp function
MINEXP(::Type{Float64}) = -7.451332191019412076235e2 # log 2^-1075
MINEXP(::Type{Float32}) = -103.97207708f0            # log 2^-150
MINEXP(::Type{Float16}) = Float16(-17.32868)         # log 2^-25

MINEXP2(::Type{Float64}) = -1075        # log 2^-1075
MINEXP2(::Type{Float32}) = -150f0       # log 2^-150
MINEXP2(::Type{Float16}) = Float16(-25) # log 2^-25

MINEXP10(::Type{Float64}) = -3.23607245338779784854769e2 # log10 2^-1075
MINEXP10(::Type{Float32}) = -45.15449934959718f0         # log10 2^-150
MINEXP10(::Type{Float16}) = Float16(-45.16)              # log10 2^-25