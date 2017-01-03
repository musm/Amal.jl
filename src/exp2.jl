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
@inline exp2_kernel{T<:LargeFloat}(x::T) = @horner_oftype(x, 1.66666666666666019037e-1,
        -2.77777777770155933842e-3,
        6.61375632143793436117e-5,
        -1.65339022054652515390e-6,
        4.13813679705723846039e-8)

# coefficients from: lib/msun/src/e_expf.c
@inline exp2_kernel{T<:SmallFloat}(x::T) = @horner_oftype(x, 1.6666625440e-1, -2.7667332906e-3)

# for values smaller than this threshold just use Taylor expansion of 1 + x*log(2)
exp2_small_thres(::Type{Float64}) = 2.0^-29
exp2_small_thres(::Type{Float32}) = 2.0f0^-14

"""
    exp2(x)

Compute the base `2` exponential of `x`, in other words ``2^x``.
"""
function exp2{T<:IEEEFloat}(x::T)
    xu = reinterpret(Unsigned, x)
    xs = xu & ~sign_mask(T)
    xsb = Int(xu >> Unsigned(8*sizeof(T)-1))

    # filter out non-finite arguments
    if xs > reinterpret(Unsigned, MAXEXP2(T))
        if xs >= exponent_mask(T)
            xs & significand_mask(T) != 0 && return T(NaN)
            return xsb == 0 ? T(Inf) : T(0.0) # exp(+-Inf)
        end
        x > MAXEXP2(T) && return T(Inf)
        x < MINEXP2(T) && return T(0.0)
    end

    # argument reduction
    if xs > reinterpret(Unsigned, T(0.5))
        if xs < reinterpret(Unsigned, T(1.5))
            if xsb == 0
                t = x - T(1.0)
                hi = t*LN2U(T)
                lo = -t*LN2L(T)
                k = 1
            else
                t = x + T(1.0)
                hi = t*LN2U(T)
                lo = -t*LN2L(T)
                k = -1
            end
        else
            n = round(x)
            k = unsafe_trunc(n)
            t = x - n
            hi = t*LN2U(T)
            lo = -t*LN2L(T)
        end
        r = hi - lo
    elseif xs < reinterpret(Unsigned, exp2_small_thres(T))
        return T(1.0) + x*T(LN2)
    else # here k = 0, no need for hi and lo, so compute approximation directly
        x *= T(LN2)
        z = x*x
        p = x - z*exp2_kernel(z)
        return T(1.0) - ((x*p)/(p - T(2.0)) - x)
    end

    # compute approximation
    z = r*r
    p = r - z*exp2_kernel(z)
    y = T(1.0) - ((lo - (r*p)/(T(2.0) - p)) - hi)
    if k > -significand_bits(T)
        # multiply by 2.0 first to prevent overflow, extending the range
        k == exponent_max(T) && return y*T(2.0)*T(2.0)^(exponent_max(T)-1)
        twopk = reinterpret(T, ((exponent_bias(T) + k) % typeof(xu)) << significand_bits(T))
        return y*twopk
    else
        # add significand_bits(T) + 1 to lift the range outside the subnormals
        twopk = reinterpret(T,
            ((exponent_bias(T) + significand_bits(T) + 1 + k) % typeof(xu)) << significand_bits(T))
        return y*twopk*T(2.0)^(-significand_bits(T) - 1)
    end
end
