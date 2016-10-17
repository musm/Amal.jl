"""
    log(x)

Compute the natural logarithm of `x`. Throws `DomainError` for negative Real
arguments. Use complex negative arguments to obtain complex results.
"""
function log end

# Based on FreeBSD lib/msun/src/e_log.c
# ====================================================
# Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved. Permission
# to use, copy, modify, and distribute this software is freely granted,
# provided that this notice is preserved.
# ====================================================

# Method :
# 1. Argument Reduction: find k and f such that
#       x = 2^k * (1+f),
#    where  sqrt(2)/2 < 1+f < sqrt(2).
#
# 2. Approximation of log(1+f). Let s = f/(2+f); based on log(1+f) = log(1+s) - log(1-s)
#       = 2s + 2/3 s**3 + 2/5 s**5 + .....,
#       = 2s + s*p
#
#    We use a special Remez algorithm on [0, f/(f+2)] (where f = sqrt(2)-1S) to
#    generate a polynomial of degree 14 to approximate p.  The maximum error of
#    this polynomial approximation is bounded by 2**-58.45. In other words,
#       p(z) = L1*s^2 + L2*s^4 + L3*s^6 + L4*s^8 + L5*s^10 + L6*s^12 + L7*s^14
#    Note that 2s = f - s*f = f - hfsq + s*hfsq, where hfsq = f*f/2.
#
#    In order to guarantee error in log below 1ulp, we compute log by
#       log(1+f) = f - s*(f - p)        (if f is not too large)
#       log(1+f) = f - (hfsq - s*(hfsq + p)). (better accuracy)
#
# 3. Finally, log(x) = k*ln2 + log(1+f).
#       = k*ln2u + (f - (hfsq - (s*(hfsq + p) + k * ln2l)))
#    Here ln2 is split into two floating point number: ln2u (upper) and ln2l (lower)


# split polynomial evaluation scheme
@inline function _log{T}(x2::T)
    x4 = x2*x2
    p1 = @horner_oftype(x4, 6.666666666666735130e-1, 2.857142874366239149e-1,
        1.818357216161805012e-1, 1.479819860511658591e-1)
    p2 = @horner_oftype(x4, 3.999999999940941908e-1, 2.222219843214978396e-1,
        1.531383769920937332e-1)
    return x2*p1 + x4*p2
end

@inline function _log{T<:SmallFloat}(x2::T)
    x4 = x2*x2
    p1 = @horner_oftype(x4, 0.6666666269302368, 0.2849878668785095)
    p2 = @horner_oftype(x4, 0.40000972151756287, 0.24279078841209412)
    return x2*p1 + x4*p2
end

function log{T}(x::T)
    x < 0 && throw(DomainError())

    # reduction, same as frexp but with custom error handling
    xu = reinterpret(Unsigned,x)
    xe = xu & exponent_mask(T)
    k = Int(xe >> significand_bits(T))

    if xe == 0 # x is subnormal
        x == 0 && return T(-Inf) # x equals +-0 
        xs = xu & sign_mask(T)
        xu $= xs
        m = leading_zeros(xu) - exponent_bits(T)
        xu <<= m
        xu $= xs
        k = 1 - m
    elseif xe == exponent_mask(T) # NaN or Inf
        return x
    end

    k -= (exponent_bias(T) - 1)
    xu = (xu & ~exponent_mask(T)) | exponent_half(T)
    f = reinterpret(T,xu)
    # in other words the above are the same as:
    # f, k = _frexp(x) # also include exception handling

    # scale up if smaller than sqrt(2)/2 for bettery accuracy
    if f < T(SQRT22)
        f *= 2
        k -= 1 
    end
    f -= 1 

    # compute approximation
    s = f/(f + 2)
    s2 = s*s
    p = _log(s2)
    hf2 = T(0.5)*f*f
    return k*LN2U(T) - ((hf2 - (s*(hf2 + p) + k*LN2L(T))) - f)
end