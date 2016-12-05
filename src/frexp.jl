"""
    frexp(val::AbstractFloat)

Return `(x, exp)` such that `x` has a magnitude in the interval ``[1/2, 1)`` or 0, and ``val = x
\\times 2^{exp}``. In other words, this return the binary significand of `val`, which when
multiplied by 2 raised to the power `exp` returns `val`.
"""
function frexp{T<:AbstractFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x, 0 # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if xs <= (~exponent_mask(T) & ~sign_mask(T)) # x is subnormal
        xs == 0 && return x, 0 # +-0
        m = leading_zeros(xs) - exponent_bits(T)
        xs <<= unsigned(m)
        xu = xs | (xu & sign_mask(T))
        k = 1 - m
    end
    k -= (exponent_bias(T) - 1)
    xu = (xu & ~exponent_mask(T)) | exponent_half(T)
    return reinterpret(T, xu), k
end
