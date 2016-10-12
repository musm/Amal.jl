tol = 5

@testset "Accuracy (max error in ulp) for $T" for T in (Float64, Float32 )
    println("Accuracy tests for $T")
    
    xx = map(T, vcat(-10:0.0002:10, -1000:0.001:1000, -120:0.023:1000, -1000:0.02:2000))
    fun_table = Dict(Amal.exp => Base.exp)
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
    fun_table = Dict(Amal.exp2 => Base.exp2)
    test_acc(T, fun_table, xx, tol)


    xx = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
    fun_table = Dict(Amal.exp10 => Base.exp10)
    test_acc(T, fun_table, xx, tol)


    # # xx = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
    # # xx = map(T, vcat(-10:0.0002:10,-10:0.0001:10, -2:0.00001:2))
    # xx = map(T, vcat(-1:0.00001:1, -1:0.000001:1,-1:0.000021:1))
    # fun_table = Dict(Amal.atan => Base.atan)
    # tol = 3
    # test_acc(T, fun_table, xx, tol)

end
