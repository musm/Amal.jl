"""
    exp(x)

Compute the natural base exponential of `x`, in other words ``e^x``.
"""
function exp end

# Method
# 1. Argument reduction: Reduce x to an r so that |r| <= 0.5*ln(2). Given x,
#    find r and integer k such that
#               x = k*ln(2) + r,  |r| <= 0.5*ln(2).
#    Here r is represented as r = hi-lo for better accuracy.
#    
# 2. Approximate exp(r) by a special rational function on [0, 0.5*ln(2)]:
#         R(r^2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r^4/360 + ...
#    
#    A special Remez algorithm on [0, 0.5*ln(2)] is used to generate a
#    polynomial to approximate R.
#        
#    The computation of exp(r) thus becomes
#                               2*r
#               exp(r) = 1 + ----------
#                             R(r) - r
#                                  r*c(r)
#                      = 1 + r + ----------- (for better accuracy)
#                                 2 - c(r)
#    where
#               c(r) = r - (P1*r^2  + P2*r^4  + ... + P5*r^10 + ...).
#     
# 3. Scale back: exp(x) = 2^k * exp(r)

# coefficients from:
# origin: FreeBSD /usr/src/lib/msun/src/e_exp.c */
# ====================================================
# Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
# Permission to use, copy, modify, and distribute this
# software is freely granted, provided that this notice
# is preserved.
# ====================================================
@inline @oftype_float function _exp{T}(hi::T, lo::T)
    r = hi - lo
    z = r*r
    p = r - z * 
    @horner(z, 1.66666666666666019037e-01,
    -2.77777777770155933842e-03,
    6.61375632143793436117e-05,
    -1.65339022054652515390e-06,
    4.13813679705723846039e-08)
    return 1.0 - ((lo - r*p/(2.0 - p)) - hi)
end

# custom coefficients
@inline @oftype_float function _exp{T<:SmallFloat}(hi::T, lo::T)
    r = hi - lo
    z = r*r
    p = r - z * 
    @horner(z, 0.1666666567325592041015625,
    -2.777527086436748504638671875e-3,
    6.451140507124364376068115234375e-5)
    return 1.0 - ((lo - r*p/(2.0 - p)) - hi)
end

@oftype_float function exp{T}(x::T)
    x > MAXEXP(T) && return Inf
    x < MINEXP(T) && return 0.0
 
    # reduce: computed as r = hi - lo for extra precision
    k = round(T(LOG2E)*x) 
    n = _trunc(k)
    hi = muladd(k, -LN2U(T), x)
    lo = k*LN2L(T)

    # compute approximation
    x = _exp(hi,lo)
    return _ldexp(x,n)
end
