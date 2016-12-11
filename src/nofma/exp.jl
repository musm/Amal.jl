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
#    Here r is represented as r = hi-lo for better accuracy.
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
# 3. Scale back: exp(x) = 2^k * exp(r)

# coefficients from: lib/msun/src/e_exp.c
@inline exp_kernel{T<:LargeFloat}(x::T) = @horner_oftype(x, 1.66666666666666019037e-1,
        -2.77777777770155933842e-3,
        6.61375632143793436117e-5,
        -1.65339022054652515390e-6,
        4.13813679705723846039e-8)

# custom coefficients
@inline exp_kernel{T<:SmallFloat}(x::T) = @horner_oftype(x, 0.1666666567325592041015625,
        -2.777527086436748504638671875e-3,
        6.451140507124364376068115234375e-5)

"""
    exp(x)

Compute the natural base exponential of `x`, in other words ``e^x``.
"""
function exp{T<:IEEEFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xsb = xu & sign_mask(T)

    # filter out non-finite arguments
    if xs > reinterpret(Unsigned, MAXEXP(T))
        if xs >= exponent_mask(T)
            if xs & significand_mask(T) != 0
                return T(NaN) 
            end
            return xsb == 0 ? T(Inf) : T(0.0) # exp(+-Inf)
        end
        x > MAXEXP(T) && return T(Inf)
        x < MINEXP(T) && return T(0.0)
    end
 
   # this implementation gives 2.7182818284590455 for exp(1.0),
   # which is well within the allowable error. however,
   # 2.718281828459045 is closer to the true value so we prefer that
   # answer, given that 1.0 is such an important argument value.
   # This is only the case when T == Float64
    if x == T(1.0) && T == Float64
        return 2.718281828459045235360
    end

    # reduce: computed as r = hi - lo for extra precision
    k = round(T(LOG2E)*x) 
    n = _trunc(k)
    hi = muladd(k, -LN2U(T), x)
    lo = k*LN2L(T)

    # compute approximation
    r = hi - lo
    z = r*r
    p = r - z * exp_kernel(z)
    x = T(1.0) - ((lo - r*p/(T(2.0) - p)) - hi)
    return _ldexp(x, n)
end
