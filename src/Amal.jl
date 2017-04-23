module Amal

# export 
#    exp, exp2, exp10,
#    log,
#    frexp, ldexp

using Base:
    significand_mask, exponent_mask, exponent_half, leading_zeros,
    Math.significand_bits, Math.exponent_bias, Math.exponent_bits,
    sign_mask, unsafe_trunc, Math.@horner, @pure

const IEEEFloat = Union{Float16,Float32,Float64}

if VERSION >= v"0.6-"
    using Base: Math.exponent_max, Math.exponent_raw_max, fpinttype
else
    # integer size of float
    fpinttype(::Type{Float64}) = UInt64
    fpinttype(::Type{Float32}) = UInt32
    fpinttype(::Type{Float16}) = UInt16

    # maximum float exponent
    @pure exponent_max{T<:IEEEFloat}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T)) - exponent_bias(T)
    # maximum float exponent without bias
    @pure exponent_raw_max{T<:IEEEFloat}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))
end


# helper functions and macros

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
include("floatmanip.jl")

include("exp.jl")
include("exp2.jl")
include("exp10.jl")

include("log.jl")

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
