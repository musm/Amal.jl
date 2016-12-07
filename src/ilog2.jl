 """
    ilog2(x)

Returns the integral part of the logarithm of `abs(x)`, using base 2 for the
logarithm. In other words, this computes the binary exponent of `x` such that
``x = significand \\times 2^exponent,`` where ``significand∈\\in [1, 2)``.

Exceptional cases (where `Int` is the machine wordsize)

* `x = 0`    returns `typemin(Int)`
* `x = ±Inf` returns `typemax(Int)`
* `x = NaN`  returns `typemax(Int)`
"""
function ilog2{T<:AbstractFloat}(x::T)
    xs = reinterpret(Unsigned, x) & ~sign_mask(T)
    xs >= exponent_mask(T) && return typemax(Int) # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && return typemin(Int)
        m = leading_zeros(xs) - exponent_bits(T)
        k = 1 - m
    end
    return k - exponent_bias(T)
end
