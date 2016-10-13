"""
    exp10(x)

Compute the base ``10`` exponential of ``x``, in other words ``10^x``.
"""
function exp10 end

# Method: see exp.jl, same idea except we do not have the hi, lo argument split

@inline @oftype_float function _exp10{T}(r::T)
    z = r*r
    p = @horner(z, 0.86858896380650363333586483349790796637535095214844,
    0.38376418216572083519366742621059529483318328857422,
    -3.3911309900330935396262077574647264555096626281738e-2,
    4.2808248822236950534292354575427452800795435905457e-3,
    -5.6858614497720487684223611424272348813246935606003e-4,
    1.74025716339103691089953973580861656955676153302193e-4,
    -4.6715889216844370609993397636117151705548167228699e-3,
    0.11757973510431163344236438206280581653118133544922,
    -1.22119911629517696738389531674329191446304321289062)
    return 1.0 + 2.0*r/(p - r)
end

@inline @oftype_float function _exp10{T<:SmallFloatTypes}(r::T)
    z = r*r
    p = @horner(z, 0.868588924407958984375,
    0.38381290435791015625,
    -5.1120765507221221923828125e-2,
    2.609136104583740234375,
    -191.65277099609375,
    6768.2021484375,
    -91928.921875)
    return 1.0 + 2.0*r/(p - r)
end

@oftype_float function exp10{T}(x::T)
    x > MAXEXP10(T) && return Inf
    x < MINEXP10(T) && return 0.0
 
    # reduce
    k = round(T(LOG210)*x) 
    n = _trunc(k)
    r = muladd(k, -LOG102U(T), x)
    r = muladd(k, -LOG102L(T), r)

    # compute approximation
    x = _exp10(r)
    return _ldexp(x,n)
end