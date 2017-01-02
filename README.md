<div align="center"><img src="https://cloud.githubusercontent.com/assets/4319522/19199051/1ecffc38-8c90-11e6-8617-19208b61a07b.jpg" alt="Amal"></img> </div>


Amal, a pure Julia math library *(work in progress)*

[![Travis Build Status](https://travis-ci.org/JuliaMath/Amal.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Amal.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/a0n8hfkv9fi59ukc/branch/master?svg=true)](https://ci.appveyor.com/project/musm/Amal-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Amal.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Amal.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Amal.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Amal.jl?branch=master)


Amal is an amalgamation of ideas from several open source math libraries, including SLEEF, Cephes, and Musl, and other references, which have all been taken into consideration for the design of this library, combining the best of breed ideas.

The Amal library principles include: avoid expensive branches, avoid table look ups, and to use FMA instructions for CPU's that support it. For CPU's with FMA instruction we take advantage of FMA and prefer polynomial functions to maximize performance if it does not sacrifice accuracy.


## Installation


We recommend running julia with `-O3` for maximal performance using `Amal.jl` and to also build a custom system image by running
```julia
julia> is_windows() && (Pkg.add("WinRPM"); using WinRPM; WinRPM.install("gcc"))
julia> include(joinpath(dirname(JULIA_HOME),"share","julia","build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then to restart `julia`; this will ensure you are taking full advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set)  if your CPU supports it.

## Usage


The exported functions presently include
```julia
exp, exp2, exp10,
log,
frexp, ldexp, ilog2
```
More functions to come in the near future.



To use  `Amal.jl`
```julia
julia> Pkg.clone("https://github.com/JuliaMath/Amal.jl.git")
julia> using Amal
julia> Amal.exp(2.0)
7.38905609893065
```