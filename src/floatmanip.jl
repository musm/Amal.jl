"""
    ldexp(x, n)

Computes ``x \\times 2^n``, using floating point manipulation, i.e. without
doing the explicit multiplitcation.

```jldoctest
julia> ldexp(5., 2)
20.0
```
"""
function ldexp{T<:IEEEFloat}(x::T, e::Integer)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && return x # +-0
        m = leading_zeros(xs) - exponent_bits(T)
        ys = xs << unsigned(m)
        xu = ys | (xu & sign_mask(T))
        k = 1 - m
        # underflow, otherwise may have integer underflow in the following n + k
        e < -50000 && return flipsign(T(0.0), x)
    end
    # For cases where e of an Integer larger than Int make sure we properly
    # overflow/underflow; this is optimized away otherwise.
    if e > typemax(Int)
        return flipsign(T(Inf), x)
    elseif e < typemin(Int)
        return flipsign(T(0.0), x)
    end
    n = e % Int
    k += n
    # overflow, if k is larger than maximum posible exponent
    if k >= exponent_raw_max(T)
        return flipsign(T(Inf), x)
    end
    if k > 0 # normal case
        xu = (xu & ~exponent_mask(T)) | (rem(k, fpinttype(T)) << significand_bits(T))
        return reinterpret(T, xu)
    else # subnormal case
        if k <= -significand_bits(T) # underflow
            # overflow, for the case of integer overflow in n + k
            e > 50000 && return flipsign(T(Inf), x)
            return flipsign(T(0.0), x)
        end
        k += significand_bits(T)
        z = T(2.0)^-significand_bits(T)
        xu = (xu & ~exponent_mask(T)) | (rem(k, fpinttype(T)) << significand_bits(T))
        return z*reinterpret(T, xu)
    end
end

"""
    frexp(val)

Return `(x, exp)` such that `x` has a magnitude in the interval ``[1/2, 1)`` or 0, and ``val = x
\\times 2^{exp}``. In other words, this return the binary significand of `val`, which when
multiplied by 2 raised to the power `exp` returns `val`.
"""
function frexp{T<:IEEEFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x, 0 # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
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

"""
    exponent(x) -> Int

Get the exponent of a normalized floating-point number. Returns the integral
part of the logarithm of `abs(x)`, using base 2 for the logarithm. In other
words, this computes the binary exponent of `x` such that ``x = significand
\\times 2^exponent,`` where ``significand \\in [1, 2)``.

Exceptional cases

* `x = 0`    returns `DomainError()`
* `x = ±Inf` returns `DomainError()`
* `x = NaN`  returns `DomainError()`
"""
function exponent{T<:IEEEFloat}(x::T)
    xs = reinterpret(Unsigned, x) & ~sign_mask(T)
    xs >= exponent_mask(T) && return throw(DomainError()) # NaN or Inf
    k = Int(xs >> significand_bits(T))
    if k == 0 # x is subnormal
        xs == 0 && throw(DomainError())
        m = leading_zeros(xs) - exponent_bits(T)
        k = 1 - m
    end
    return k - exponent_bias(T)
end

"""
    significand(x)

Extract the `significand(s)` (a.k.a. mantissa), in binary representation, of a
floating-point number. If `x` is a non-zero finite number, then the result will be
a number of the same type on the interval ``[1,2)``. Otherwise `x` is returned.

```jldoctest
julia> significand(15.2)/15.2
0.125

julia> significand(15.2)*8
15.2
"""
function significand{T<:IEEEFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xs >= exponent_mask(T) && return x # NaN or Inf
    if xs <= (~exponent_mask(T) & ~sign_mask(T)) # x is subnormal
        xs == 0 && return x # +-0
        m = unsigned(leading_zeros(xs) - exponent_bits(T))
        xs <<= m
        xu = xs | (xu & sign_mask(T))
    end
    xu = (xu & ~exponent_mask(T)) | exponent_one(T)
    return reinterpret(T, xu)
end
