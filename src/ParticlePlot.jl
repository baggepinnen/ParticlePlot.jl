module ParticlePlot

export kde, heatBoxPlot, densityplot, pplot, pploti
using RecipesBase, StatsBase
"""
Weighted kernel density estimate of the data `x` ∈ ℜN with weights `w` ∈ ℜN
`xi, densityw, density = kde(x,w)`
The number of grid points is chosen automatically and will approximately be equal to N/3
The bandwidth of the Gaussian kernel is chosen based on Silverman's rule of thumb
returns both weighted and non-weighted densities `densityw, density`
"""
function kde(x,w)
    @assert all(w .>= 0) "All weights must be non-negative"
    x   = x[:]
    e   = StatsBase.histrange(x,ceil(Int,length(x)/3))
    nb  = length(e)-1
    np  = length(x)
    I   = sortperm(x)
    s   = (rand()/nb+0):1/nb:e[end]
    j   = zeros(Float64,nb)
    bo = 1
    ilast = 0
    for i = 1:np
        ii = I[i]
        for b = bo:nb
            if x[ii] <= e[b+1]
                j[b] += w[ii]
                bo = b
                ilast = i
                break
            end
        end
    end
    j[end] += sum(w[I[ilast+1:np]])
    @assert all(j .>= 0) "j"
    xi      = e[1:end-1]+0.5(e[2]-e[1])
    σ       = std(x)
    if σ == 0 # All mass on one point
        return x,ones(x),ones(x)
    end
    h        = 1.06σ*np^(-1/5)
    K(x)     = exp(-(x/h).^2/2) / √(2π)
    densityw = [1/h*sum(j.*K(xi.-xi[i])) for i=1:nb] # This density is effectively normalized by nb due to sum(j) = 1
    density  = [1/(np*h)*sum(K(x.-xi[i])) for i=1:nb]

    # @assert all(densityw .>= 0) "densityw"
    # @assert all(density  .>= 0) "density"

    return xi, densityw, density

end

@userplot DensityPlot
@recipe function f(dp::DensityPlot)
    x = dp.args[1]
    w = length(dp.args) > 1 ? dp.args[2] : ones(x)/length(x)
    xi, densityw, density = kde(x,w)
    if maximum(x)-minimum(x) > 0
        title --> "Kernel density estimate"
        seriestype := :path
        @series begin
            label --> "Weighted density"
            linecolor --> :blue
            xi,densityw
        end
        @series begin
            label --> "Non-weighted density"
            linecolor --> :red
            xi,density
        end
    end
end

function heatBoxPlot(plt, x, t, minmax, nbinsy=30)
    if maximum(x)-minimum(x) > 0
        heatmap!(plt, x=(t-1):1/length(x):t,y=x',nbins=(1,nbinsy), ylims=tuple(minmax...))
    end
    return
end


""" `pplot(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata)`
To be called inside a particle filter, plots either particle density (`density=true`) or individual particles (`density=false`) \n
Will plot all the real states in `xIndices` as well as the expected vs real measurements of `yIndices`.
Arguments: \n
* `x`: `Array(M,N)`. The states for each patyicle where `M` number of states, `N` number of Particles
* `w`: `Array(N)`. weight of each particle
* `y`: `Array(R,T)`. All true outputs. `R` is number of outputs, `T` is total number of time steps (will only use index `t`)
* `yhat`: `Array(R,N)` The expected output per particle. `R` is number of outputs, `N` number of Particles
* `N`, Number of particles
* `a`, `Array(N)`, reorderng of particles (e.g. `1:N`)
* `t`, Current time step
* `xreal`: `Array(M,T)`. All true states. `R` is number of states, `T` is total number of time steps (will only use index `t`)
* `xhat`: Not used
* `xOld`: Same as `x`, but for previous time step, only used when `!density` to show states origins
* `pdata`: Persistant data for plotting. Set to void in first call and pdataOut on remaining \n
* `density = true` To only plot the particle trajectories, set (`leftOnly=false`)\n
* `leftOnly = false`\n
* `xIndices = 1:size(x,1)`\n
* `yIndices = 1:size(y,1)`\n
Returns: `pdataOut`
"""
function pplot(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata; density = true, leftOnly = false, xIndices = 1:size(x,1), yIndices = 1:size(y,1), slidef=0.9)
    immerse()
    cols = leftOnly?1:2
    grd = (r,c) -> (r-1)*cols+c
    println("Surviving: "*string((N-length(setdiff(Set(1:N),Set(a))))/N))
    vals = [x;yhat]
    realVals = [xreal;y]
    if !density
        valsOld = [xOld;yhat]
    end

    pltIdx = [xIndices; size(x,1)+yIndices]
    if pdata == Void
        pdata = (subplot(layout=cols*ones(Int,length(pltIdx))), zeros(length(pltIdx),2))
        gui(pdata[1])
    end
    p, minmax = pdata
    dataMin = minimum(vals[pltIdx,:],2)
    dataMax = maximum(vals[pltIdx,:],2)
    minmax = [min(minmax[:,1], dataMin)*slidef+(1-slidef)*dataMin max(minmax[:,2], dataMax)*slidef+(1-slidef)*dataMax]
    #c = (w[:]-minimum(w))*3
    c = w[:]*5*N

    for (i, val) in enumerate(pltIdx)
        if !leftOnly
            oldFig = p.plts[grd(i,2)].o[1]
            newPlot = plot()
        end
        #Plot the heatmap on the left plot
        density && heatBoxPlot(p.plts[grd(i,1)], vals[val,:], t, minmax[i,:])

        for j = 1:N
            if !density
                #Plot the line on the left plot
                plot!(p.plts[grd(i,1)], [t-1.5,t-1], [valsOld[val,j], valsOld[val,j]], legend=false)
                plot!(p.plts[grd(i,1)], [t-1,t-0.5], [valsOld[val,a[j]], vals[val,j]], legend=false)
            end
            #Plot each of the dots on the right side
            !leftOnly && !density && plot!(newPlot, [j,j], [vals[val,j]',vals[val,j]'], marker=(:circle, :red, c[j])  , legend=false)
        end
        if !leftOnly
            density && densityplot!(newPlot, vals[val,:], w , ylims = tuple(minmax[i,:]...), c=:blue)
            #Fix the bug with updating plot
            p.plts[grd(i,2)] = copy(newPlot)
            p.plts[grd(i,2)].o = (oldFig, p.plts[grd(i,2)].o[2])
        end
        #Plot Real State Here
        plot!(p.plts[grd(i,1)], (t-1):t, [realVals[val,t], realVals[val,t]], legend=false, color=:black, linewidth=5)
    end
    gui(p)
    (p, minmax)
end

# """ `pploti(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata)`
# Same function as pplot but with options to skip, wait and quit.
# """
function pploti(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata; density = true, leftOnly = false, xIndices = 1:size(x,1), yIndices = 1:size(y,1), slidef=0.9)
    cols = leftOnly?1:2
    pltIdx = [xIndices; size(x,1)+yIndices]
    if pdata == Void
        pdata = (subplot(layout=cols*ones(Int,length(pltIdx))), Array{Float64,2}(length(pltIdx),2),0,0)
        gui(pdata[1])
    end
    p, minmax, skip, skipped = pdata
    if skip > 0 && skipped < skip
        skipped += 1
        return (p, minmax, skip, skipped)
    elseif skip > 0
        skipped = 0
        skip = 0
    end
    p, minmax = pplot(x, w, y, yhat, N, a, t, xreal, xhat, xOld, (p, minmax), density = density, leftOnly = leftOnly, xIndices = xIndices, yIndices = yIndices)
    print("Waiting for command. q to Quit, ^D to run all, s NN to skip NN steps:\n")
    line = readline(STDIN)
    if line == "q\n"
        error("Quitting")
    elseif contains(line, "s")
        ss = split(strip(line,'\n'))
        skip = parse(Int,ss[2])
    end
    (p, minmax, skip, skipped)
end
end # module
