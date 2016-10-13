module Amal

using Base: significand_bits, exponent_bias, exponent_mask, exponent_one, Math.@horner

typealias FloatTypes Union{Float16,Float32,Float64}
typealias IntTypes Union{Int16,Int32,Int64}

typealias LargeFloatTypes Union{Float64}
typealias SmallFloatTypes Union{Float16,Float32}

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

floattype(::Type{Int64}) = Float64
floattype(::Type{Int32}) = Float32
floattype(::Type{Int16}) = Float16

# convert float to corresponding integer type of the same size
asint{T<:FloatTypes}(::Type{T}, x) = inttype(T)(x)

# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IntTypes}(m::T) = reinterpret(floattype(T), m << significand_bits(floattype(T)))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:FloatTypes}(d::T) = reinterpret(inttype(T), d) >> significand_bits(T)

exponent_max{T<:FloatTypes}(::Type{T}) = asint(T, exponent_mask(T) >> significand_bits(T))

# import Base: unsafe_trunc
# unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x)) # Float16 hack for v0.5
# unsafe truncate x and return integer of the same size as x
_trunc{T<:FloatTypes}(x::T) = unsafe_trunc(inttype(T), x)

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
