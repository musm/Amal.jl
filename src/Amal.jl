module Amal

# export 
#    exp, exp2, exp10,
#    log,
#    frexp, ldexp

using Base:
    significand_mask, Math.significand_bits,
    exponent_mask, exponent_half,
    Math.exponent_bias, Math.exponent_bits, Math.exponent_max,
    sign_mask, unsafe_trunc, leading_zeros, fpinttype, Math.@horner

const IEEEFloat = Union{Float16,Float32,Float64}

# helper functions and macros

# evaluate p[1] + x * (p[2] + x * (....)), i.e. a polynomial via Horner's rule
# and convert coefficients to same type as x
macro horner_oftype(x, p...)
    ex = :(oftype($(esc(x)),$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd(t, $ex, oftype($(esc(x)),$(esc(p[i])))))
    end
   return Expr(:block, :(t = $(esc(x))), ex)
end

#####################################################

include("constants.jl")

include("floatmanip.jl")

include("exp.jl")
include("exp10.jl")

include("exp2.jl")
include("log.jl")

# Float16 definitions
for func in (:exp, :exp2, :exp10, :log)
    @eval begin
        $func(a::Float16) = Float16($func(Float32(a)))
    end
end
ldexp(x::Float16, q::Integer) = Float16(ldexp(Float32(x), q))

# fallback definitions
for func in (:exp, :exp2, :exp10, :log)
    @eval begin
        $func(x::Real) = $func(float(x))
    end
end
end
