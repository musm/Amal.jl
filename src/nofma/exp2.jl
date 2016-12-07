"""
    exp2(x)

Compute the base `2` exponential of `x`, in other words ``2^x``.
"""
function exp2 end

# Based on FreeBSD lib/msun/src/e_exp.c
# ====================================================
# Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved. Permission
# to use, copy, modify, and distribute this software is freely granted,
# provided that this notice is preserved.
# ====================================================

# Method
# 1. Argument reduction: Reduce x to an r so that |r| <= 0.5*ln(2). Given x,
#    find r and integer k such that
#       x = k*ln(2) + r,  |r| <= 0.5*ln(2).
#    Here r is represented as r = hi - lo for better accuracy.
#    
# 2. Approximate exp(r) by a special rational function on [0, 0.5*ln(2)]:
#       R(r^2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r^4/360 + ...
#    
#    A special Remez algorithm on [0, 0.5*ln(2)] is used to generate a
#    polynomial to approximate R.
#        
#    The computation of exp(r) thus becomes
#                       2*r
#       exp(r) = 1 + ----------
#                     R(r) - r
#                          r*c(r)
#              = 1 + r + ----------- (for better accuracy)
#                         2 - c(r)
#    where
#       c(r) = r - (P1*r^2  + P2*r^4  + ... + P5*r^10 + ...).
#      
# 4. To obtain exp2 we simply scale the input argument by log(2) and use
#    the coefficients from exp(x) Note: trying the same for exp10, results
#    in a very innacurate function

# coefficients from: lib/msun/src/e_exp.c
@inline _exp2{T}(x::T) = @horner_oftype(x, 1.66666666666666019037e-1,
        -2.77777777770155933842e-3,
        6.61375632143793436117e-5,
        -1.65339022054652515390e-6,
        4.13813679705723846039e-8)

# custom coefficients
@inline _exp2{T<:SmallFloat}(x::T) = @horner_oftype(x, 0.1666666567325592041015625,
        -2.777527086436748504638671875e-3,
        6.451140507124364376068115234375e-5)

function exp2{T}(x::T)
    x > MAXEXP2(T) && return T(Inf)
    x < MINEXP2(T) && return T(0.0)
 
    # reduce: computed as r = hi - lo for extra precision
    k = round(x)
    n = _trunc(k)
    t = x - k
    hi = t*LN2U(T)
    lo = -t*LN2L(T)

    # compute approximation
    r = hi - lo
    z = r*r
    p = r - z * _exp2(z)
    x = T(1.0) - ((lo - r*p/(T(2.0) - p)) - hi)
    return _ldexp(x, n)
end
