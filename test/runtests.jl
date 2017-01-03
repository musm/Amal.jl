using Amal
using Suppressor
using Base.Test

function countulp{T}(::Type{T}, X::AbstractFloat,Y::AbstractFloat)
    x,y = T(X),T(Y)
    (isnan(x) && isnan(y)) && return zero(T)
    (isnan(x) || isnan(y)) && return T(1000)
    (isinf(x) && isinf(y)) && return (signbit(x) == signbit(y)) ? zero(T) : T(1001)
    (x == zero(typeof(X)) && y == zero(typeof(Y))) && return (signbit(x) == signbit(y)) ? zero(T) : T(1002)
    if isfinite(x) && isfinite(y)
        d = abs(X - Y)
        return T(d/eps(y))
        # alternative definition
        # k = abs(frexp(Y)[2])
        # return T(ldexp(d,Base.significand_bits(T)-k+1))
    end
    return T(1003)
end
countulp{T<:AbstractFloat}(x::T, y::T) = countulp(T,x,y)

function cmpdenorm{T}(::Type{T}, X::AbstractFloat, Y::AbstractFloat)
    x,y = T(X),T(Y)
    (isnan(x) && isnan(y)) && return signbit(x) == signbit(y)
    (isinf(x) && isinf(y)) && return signbit(x) == signbit(y)
    (isfinite(x) && isfinite(y)) && return signbit(x) == signbit(y)
    return false
end

function _test_acc{T}(::Type{T}, fun_test::Function, fun_ref::Function, xx)
    vmax  = zero(T)
    vmean = zero(T)
    xmax = map(zero, first(xx))

    local vtest::T
    local vtrue::BigFloat
    local u::T
    @inbounds for x in xx
        vtest = fun_test(x...)
        vtrue = fun_ref(map(BigFloat,x)...)
        u = countulp(T,vtest,vtrue)

        vmax = max(vmax, u)
        xmax = vmax == u ? x : xmax
        vmean += u
    end
    vmean = vmean/length(xx)
    vmax, vmean, xmax
end

strip_module_name(f::Function) = last(split(string(f), '.')) # strip function name from qualified name

# test the accuracy of a function where fun_table is a Dict mapping the function you want
# to test to a reference function
# xx is an array of values (which may be tuples for multiple arugment functions)
# tol is the acceptable tolerance to test against
function test_acc(T::Type, fun_table, xx, tol)
    @testset "accuracy $(strip_module_name(fun_test))" for (fun_test, fun_ref) in fun_table
        vmax, vmean, xmax = _test_acc(T, fun_test, fun_ref, xx)
        @test trunc(vmax,2) <= tol
        # print test result
        fmtxloc = isa(xmax, Tuple) ? string('(', join((@sprintf("%.5f", x) for x in xmax), ", "), ')') : @sprintf("%.5f", xmax)
        println(rpad(strip_module_name(fun_test), 18, " "), ": max ", @sprintf("%0.8f", vmax),
            rpad(" at x = "*fmtxloc, 40, " "), ": mean ", @sprintf("%0.8f", vmean))
    end
end

function runtests()
    @testset "Amal" begin
    include("exceptional.jl")
    include("accuracy.jl")
    end
end

runtests()
