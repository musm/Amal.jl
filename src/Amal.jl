module Amal

export exp, exp2, exp10,
       log

typealias IntrinsicFloats Union{Float16,Float32,Float64}
typealias IntTypes Union{Int16,Int32,Int64}

typealias LargeFloat Union{Float64}
typealias SmallFloat Union{Float16,Float32}


using Base: significand_bits, exponent_bias, exponent_mask

#################################################################################################

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

exponent_max{T<:IntrinsicFloats}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))
# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IntrinsicFloats}(::Type{T}, m::Integer) = reinterpret(T, (m % inttype(T)) << significand_bits(T))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:IntrinsicFloats}(d::T) = (reinterpret(inttype(T), d) >> significand_bits(T)) % Int

# unsafe truncate x
_trunc{T<:IntrinsicFloats}(x::T) = unsafe_trunc(Int, x)
# hack: until we have better Float16 support
_trunc(x::Float16) = (isnan(x) || isinf(x)) ? typemin(Int16) : trunc(Int, x)

#################################################################################################

include("macros.jl")
include("constants.jl")

include("ldexp.jl")

if FMA_FAST
    include("poly/exp.jl")   # slightly less accurate than rational version on non FMA systems
    include("poly/exp2.jl")  # more accurate than rational version for FMA systems
else
    include("exp.jl")
    include("exp2.jl") # more accurate for non FMA
end
# include("exp10.jl") # try to develop a better rational approximation
include("poly/exp10.jl") # better than rational approx for FMA and non FMA systems

include("log.jl")

end
