module Amal

using Base: significand_bits, exponent_bias, exponent_mask, exponent_one

typealias FloatTypes Union{Float16,Float32,Float64}
typealias IntTypes Union{Int16,Int32,Int64}

typealias LargeFloatTypes Union{Float64}
typealias SmallFloatTypes Union{Float16,Float32}

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

floattype(::Type{Int64}) = Float64
floattype(::Type{Int32}) = Float32
floattype(::Type{Int16}) = Float16

# convert float to corresponding integer type of the same size
asint{T<:FloatTypes}(::Type{T}, x) = inttype(T)(x)

# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IntTypes}(m::T) = reinterpret(floattype(T), m << significand_bits(floattype(T)))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:FloatTypes}(d::T) = reinterpret(inttype(T), d) >> significand_bits(T)

exponent_max{T<:FloatTypes}(::Type{T}) = asint(T, exponent_mask(T) >> significand_bits(T))

# import Base: unsafe_trunc
# unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x)) # Float16 hack for v0.5
# unsafe truncate x and return integer of the same size as x
_trunc{T<:FloatTypes}(x::T) = unsafe_trunc(inttype(T), x)



# evaluate p[1] + x * (p[2] + x * (....)), i.e. a polynomial via Horner's rule
# and convert coefficients to same type as x
macro horner_oftype(x, p...)
    @gensym val
    ex = :(oftype($x,$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd($val, $ex, oftype($x,$(esc(p[i])))))
    end
   return Expr(:block, :($val = $(esc(x))), ex)
end

# Similar to @horner, but split into even and odd coefficients.
macro horner_split_oftype(x,p...)
    @gensym t1
    @gensym t2
    blk = quote
        $t1 = $(esc(x))
        $t2 = $(esc(x)) * $(esc(x))
    end
    n = length(p)
    p0 = :(oftype($x,$(esc(p[1]))))
    if isodd(n)
        ex_o = :(oftype($x,$(esc(p[end-1]))))
        ex_e = :(oftype($x,$(esc(p[end]))))
        for i = n-3:-2:2
            ex_o = :(muladd($(t2), $ex_o, oftype($x,$(esc(p[i])))))
        end
        for i = n-2:-2:2
            ex_e = :(muladd($(t2), $ex_e, oftype($x,$(esc(p[i])))))
        end
    elseif iseven(n)
        ex_o = :(oftype($x,$(esc(p[end]))))
        ex_e = :(oftype($x,$(esc(p[end-1]))))
        for i = n-2:-2:2
            ex_o = :(muladd($(t2), $ex_o, oftype($x,$(esc(p[i])))))
        end
        for i = n-3:-2:2
            ex_e = :(muladd($(t2), $ex_e, oftype($x,$(esc(p[i])))))
        end
    end
    push!(blk.args,:($(p0) + $(t1)*$(ex_o) + $(t2)*$(ex_e)))
    blk
end

function _numeric(T::Symbol, ex::Union{Symbol,Expr}) # ex is either a Symbol or Expr
    isa(ex, Symbol) && return ex
    if isa(ex, Expr)
        if is(ex.head, :line)
            # skip line information
            return ex
        else
            # warning Inf is not correctly parsed as a numeric type, so we have to special case this
            ex.args = map(x -> (isa(x, Number) || (x == :Inf) || (x == :Inf32) || 
                            (x == :Inf16)) ? :($T($x)) : _numeric(T,x), ex.args)        
        end
        return ex
    end
end

# wrap signature to make the function parameteric over the float constants
# e.g. f(x::AbstractFloat) = x + 1.0  -> f{T<:AbstractFloat}(x::T) = x + T(1.0)
macro oftype_float(ex::Expr)
    ex_sig = ex.args[1]
    ex_body = ex.args[2]
    # check that ex is a function in the form  `f(x) = ... `  or `function f(x) ... end`
    if (is(ex.head, :(=)) || is(ex.head, :function)) && is(ex_sig.head, :call)
        # if signature is in the form f(x::Type)
        if isa(ex_sig.args[1], Symbol) && isa(ex_sig.args[2], Expr)
            type_parameter = :T
             # assume first type signature is the same for all other arguments
            type_signature = ex_sig.args[2].args[2]
            function_name = ex_sig.args[1]
            ex_sig.args[1] = Expr(:curly, function_name, Expr(:<:, type_parameter ,type_signature))
            for i = 2:length(ex_sig.args)  # copy new type signature to all other arguments
                ex_sig.args[i].args[2] = type_parameter
            end
            ex_body  = _numeric(type_parameter, ex_body)
            return esc(ex)
        # if signature is in the form f{T<:Type}(x::T)
        elseif is(ex_sig.args[1].head, :curly) && isa(ex_sig.args[2],Expr)
             # assume first argument type signature is the same for all other numeric type signatures arguments
            type_parameter = ex_sig.args[2].args[2]
            ex_body  = _numeric(type_parameter, ex_body)
            return esc(ex)
        # if signature is in the form f{T}(x::T)
        elseif is(ex_sig.args[1].head, :curly) 
            type_parameter = ex_sig.args[1].args[2]
            ex_body  = _numeric(type_parameter, ex_body)
            return esc(ex)
        end
    end
end

include("constants.jl")

include("ldexp.jl")

include("exp.jl") # exp, exp2
# include("poly/exp.jl") # slightly less accurate than rational approximation

# include("exp10.jl")
include("poly/exp10.jl") # slightly more accurate than rational approximation

# include("poly/exp2.jl")


end
