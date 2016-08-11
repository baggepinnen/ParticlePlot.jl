# ParticlePlot

[![Build Status](https://travis-ci.org/baggepinnen/ParticlePlot.jl.svg?branch=master)](https://travis-ci.org/baggepinnen/ParticlePlot.jl)

[![Coverage Status](https://coveralls.io/repos/baggepinnen/ParticlePlot.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/baggepinnen/ParticlePlot.jl?branch=master)

[![codecov.io](http://codecov.io/github/baggepinnen/ParticlePlot.jl/coverage.svg?branch=master)](http://codecov.io/github/baggepinnen/ParticlePlot.jl?branch=master)


`pplot(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata)`
To be called inside a particle filter, plots either particle density (`density=true`) or individual particles (`density=false`) 
Will plot all the real states in `xIndices` as well as the expected vs real measurements of `yIndices`.
Arguments:
* `x`: `Array(M,N)`. The states for each particle where `M` number of states, `N` number of Particles
* `w`: `Array(N)`. weight of each particle
* `y`: `Array(R,T)`. All true outputs. `R` is number of outputs, `T` is total number of time steps (will only use index `t`)
* `yhat`: `Array(R,N)` The expected output per particle. `R` is number of outputs, `N` number of Particles
* `N`, Number of particles
* `a`, `Array(N)`, reordering of particles (e.g. `1:N`)
* `t`, Current time step
* `xreal`: `Array(M,T)`. All true states. `R` is number of states, `T` is total number of time steps (will only use index `t`)
* `xhat`: Not used
* `xOld`: Same as `x`, but for previous time step, only used when `!density` to show states origins
* `pdata`: Persistant data for plotting. Set to void in first call and pdataOut on remaining 
* `density = true` To only plot the particle trajectories, set (`leftOnly=false`)
* `leftOnly = false`\n
* `xIndices = 1:size(x,1)`
* `yIndices = 1:size(y,1)`
Returns: `pdataOut`
