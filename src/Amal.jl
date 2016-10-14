module Amal

export exp, exp2, exp10

using Base: significand_bits, exponent_bias, exponent_mask, exponent_one, Math.@horner

typealias IntrinsicFloats Union{Float16,Float32,Float64}
typealias IntTypes Union{Int16,Int32,Int64}

typealias LargeFloat Union{Float64}
typealias SmallFloat Union{Float16,Float32}

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

floattype(::Type{Int64}) = Float64
floattype(::Type{Int32}) = Float32
floattype(::Type{Int16}) = Float16

# convert float to corresponding integer type of the same size
asint{T<:IntrinsicFloats}(::Type{T}, x) = inttype(T)(x)

# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IntTypes}(m::T) = reinterpret(floattype(T), m << significand_bits(floattype(T)))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:IntrinsicFloats}(d::T) = reinterpret(inttype(T), d) >> significand_bits(T)

exponent_max{T<:IntrinsicFloats}(::Type{T}) = asint(T, exponent_mask(T) >> significand_bits(T))

# unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x)) # Float16 hack for v0.5
# unsafe truncate x and return integer of the same size as x
_trunc{T<:IntrinsicFloats}(x::T) = unsafe_trunc(inttype(T), x)
# hack: until we have native Float16 support
_trunc(x::Float16) = (isnan(x) || isinf(x)) ? typemin(Int16) : trunc(Int16, x)

include("macros.jl")
include("constants.jl")

include("ldexp.jl")

# include("exp.jl")
include("poly/exp.jl") # slightly less accurate than rational approximation on non fma systems

# include("exp10.jl") # try to develop a better rational approximation
include("poly/exp10.jl") # better than rational approx for fma and non fma

# include("exp2.jl") # more accurate for non fma
include("poly/exp2.jl")  # more accurate than rational for fma

end
