module Amal

export exp, exp2, exp10,
       log, ilog2,
       frexp

using Base: significand_mask, significand_bits, exponent_bias, exponent_mask,
    exponent_half, leading_zeros, exponent_bits, sign_mask

typealias IEEEFloat Union{Float16,Float32,Float64}
typealias LargeFloat Union{Float64}
typealias SmallFloat Union{Float16,Float32}

#####################################################
# helper functions and macros

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
_div{T<:Base.BitSigned64}(x::T, y::T) = 
    Base.llvmcall("%3 = sdiv i64 %0, %1 ret i64 %3", T, Tuple{T, T}, x, y)


function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = 
        $(muladd(nextfloat(T(1.0)), nextfloat(one(T)), -nextfloat(T(1.0), 2)) != zero(T))
end
const IS_FMA_FAST = is_fma_fast(Float64) && is_fma_fast(Float32)

# evaluate p[1] + x * (p[2] + x * (....)), i.e. a polynomial via Horner's rule
# and convert coefficients to same type as x
macro horner_oftype(x, p...)
    @gensym val
    ex = :(oftype($(esc(x)),$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd($val, $ex, oftype($(esc(x)),$(esc(p[i])))))
    end
   return Expr(:block, :($val = $(esc(x))), ex)
end

#####################################################

include("constants.jl")

# floating point manipulation functions
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

# Float16 definitions
for func in (:exp,:exp2,:exp10,:log)
    @eval begin
        $func(a::Float16) = Float16($func(Float32(a)))
    end
end
ldexp(x::Float16, q::Integer) = Float16(ldexp(Float32(x), q))

for func in (:exp,:exp2,:exp10,:log)
    @eval begin
        $func(x::Real) = $func(float(x))
    end
end
end
