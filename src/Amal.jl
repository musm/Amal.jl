module Amal

export exp, exp2, exp10,
       log, ilog2,
       frexp

typealias IEEEFloat Union{Float16,Float32,Float64}

typealias LargeFloat Union{Float64}
typealias SmallFloat Union{Float16,Float32}


using Base: significand_mask, significand_bits, exponent_bias, exponent_mask,
    exponent_half, leading_zeros, exponent_bits, sign_mask

#################################################################################################

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

exponent_raw_max{T<:IEEEFloat}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))
# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IEEEFloat}(::Type{T}, m::Integer) = reinterpret(T, (m % inttype(T)) << significand_bits(T))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:IEEEFloat}(d::T) = reinterpret(Signed, d) >> significand_bits(T)

# unsafe truncate x
_trunc{T<:IEEEFloat}(x::T) = unsafe_trunc(Int, x)
# hack: until we have better Float16 support
_trunc(x::Float16) = (isnan(x) || isinf(x)) ? typemin(Int16) : trunc(Int, x)

# unsafe div
_div{T<:Base.BitSigned64}(x::T, y::T) = Base.llvmcall("%3 = sdiv i64 %0, %1 ret i64 %3", T, Tuple{T, T}, x, y)
#################################################################################################

include("macros.jl")
include("constants.jl")

#Floating point manipulation functions
include("frexp.jl")
include("ldexp.jl")

if IS_FMA_FAST
    include("exp.jl")   # slightly less accurate than rational version on non FMA systems
    include("exp2.jl")  # more accurate than rational version for FMA systems
else
    include("nofma/exp.jl")
    include("nofma/exp2.jl") # more accurate for non FMA systems
end
include("exp10.jl")

include("log.jl")
include("ilog2.jl")

for f in (:exp,:exp2,:exp10,:log)
    @eval begin
        ($f)(x::Real) = ($f)(float(x))
    end
end

end
