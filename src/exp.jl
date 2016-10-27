"""
    exp(x)

Compute the natural base exponential of `x`, in other words ``e^x``.
"""
function exp end

#  Method
#    1. Argument reduction: Reduce x to an r so that |r| <= 0.5*ln(2). Given x,
#       find r and integer k such that
#
#                x = k*ln(2) + r,  |r| <= 0.5*ln(2).
#
#    2. Approximate exp2(r) by a polynomial on the interval [-0.5*ln(2), 0.5*ln(2)]:
#
#           exp(x) = 1.0 + x + polynomial(x),
#                  = polynomial(x) + x + 1 (for better accuracy if no fma)
#
#    3. Scale back: exp(x) = 2^k * exp(r)

@inline _exp{T}(x::T) = @horner_oftype(x, 1.0, 1.0, 0.5,
    0.16666666666666685170383743752609007060527801513672,
    4.1666666666666692109277647659837384708225727081299e-2,
    8.3333333333159547579027659480743750464171171188354e-3,
    1.38888888888693412537733706813014578074216842651367e-3,
    1.9841269898657093212653024227876130680670030415058e-4,
    2.4801587357008890921336585755341275216778740286827e-5,
    2.7557232875898009206386968239499424271343741565943e-6,
    2.7557245320026768203034231441428403286408865824342e-7,
    2.51126540120060271373185023340013355408473216812126e-8,
    2.0923712382298872819985862227861600493028504388349e-9)

@inline _exp{T<:SmallFloat}(x::T) = @horner_oftype(x, 1.0, 1.0, 0.5,
    0.1666666567325592041015625,
    4.1666455566883087158203125e-2,
    8.333526551723480224609375e-3,
    1.39357591979205608367919921875e-3,
    1.97799992747604846954345703125e-4)

function exp{T}(x::T)
    x > MAXEXP(T) && return T(Inf)
    x < MINEXP(T) && return T(0.0)
    
    # reduce
    k = round(T(LOG2E)*x)
    n = _trunc(k)
    r = muladd(k, -LN2U(T), x)
    r = muladd(k, -LN2L(T), r)

    # compute approximation
    u = _exp(r)
    return _ldexp(u,n)
end