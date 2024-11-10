function get_pDf(x0,y0,P,D,d)
    p = @. ((d*x0^2+d)*y0+D*x0^2-P*x0)/((P*x0+D)*y0^2+d*y0*x0^2+D*x0^2)
    Df = @. ((P*x0+D)*y0^2+d*y0*x0^2+D*x0^2)/((x0^2+1)*y0^2+x0^2)
    return p, Df
end

function get_fitted_PDd(x0,y0,P,D,d)
    p, Df = get_pDf(x0,y0,P,D,d)
    Pf = @. Df*x0*(1-p)
    df = @. Df*y0*p
    return Pf, Df, df
end
export get_fitted_PDd

function jacobian(pars,P,D,d)
    ns = length(P)
    x0 = pars[1]
    y0 = pars[2]
    p, Df = get_pDf(x0,y0,P,D,d)
    ddPdx0 = @. Df*(1-p)
    ddPdy0 = 0
    ddPdp = -Df.*x0
    ddPdDf = @. x0*(1-p)
    ddDdx0 = fill(0,ns)
    ddDdy0 = fill(0,ns)
    ddDdp = fill(0,ns)
    ddDdDf = fill(1,ns)
    ddddx0 = fill(0,ns)
    ddddy0 = Df.*p
    ddddp = Df.*y0
    ddddDf = y0.*p
    J = zeros(3*ns,2*(ns+1))
    iP = 1:ns
    iD = ns+1:2*ns
    id = 2*ns+1:3*ns
    jx0 = 1
    jy0 = 2
    jp = 3:2+ns
    jDf = 3+ns:2+2*ns
    J[iP,jx0] .= ddPdx0
    J[iP,jy0] .= ddPdy0
    J[iP,jp] .= diagm(ddPdp)
    J[iP,jDf] .= diagm(ddPdDf)
    J[iD,jx0] .= ddDdx0
    J[iD,jy0] .= ddDdy0
    J[iD,jp] .= diagm(ddDdp)
    J[iD,jDf] .= diagm(ddDdDf)
    J[id,jx0] .= ddddx0
    J[id,jy0] .= ddddy0
    J[id,jp] .= diagm(ddddp)
    J[id,jDf] .= diagm(ddddDf)
    return J
end

function plot(x0,y0,E,
              P::AbstractVector,
              D::AbstractVector,
              d::AbstractVector;
              xlim = [0,maximum([P./D;x0+sqrt(E[1,1])])],
              ylim = [0,maximum([d./D;y0+sqrt(E[2,2])])],
              legend = false,
              nsigma = 2,
              plot_options...)
    nstep = 50
    x = range(xlim[1],xlim[2],nstep)
    y = @. y0 - x*y0/x0
    J = zeros(nstep,2)
    J[:,1] = x.*y0/x0^2
    J[:,2] = 1 .- x/x0
    covmat = J * E * transpose(J)
    sy = sqrt.(diag(covmat))
    p = Plots.plot(x,y,ribbon=nsigma*sy;legend=legend,xlim=xlim,ylim=ylim,plot_options...)
    Plots.plot!(P./D,d./D;seriestype=:scatter,legend=legend,plot_options...)
    Plots.plot!([0,x0],[y0,0];seriescolor=:black,legend=legend)
    Plots.xlabel!("P/D")
    Plots.ylabel!("d/D")
    return p
end
export plot

function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector)
    function residuals(par)
        Pf, Df, df = get_fitted_PDd(par[1],par[2],P,D,d)
        return [Pf.-P;Df.-D;df.-d]
    end
    function misfit(par::AbstractVector)
        r = residuals(par)
        return sum(r.^2)
    end
    X = P./D
    Y = d./D
    slope = (minimum(Y)-maximum(Y))/(maximum(X)-minimum(X))
    x0i = maximum(X) - minimum(Y)/slope
    y0i = maximum(Y) - slope*minimum(X)
    fit = Optim.optimize(misfit,[x0i,y0i])
    pars = Optim.minimizer(fit)
    x0 = pars[1]
    y0 = pars[2]
    ns = length(P)
    s2 = misfit(pars)/(ns-2)
    s2 = misfit(pars)/(ns-2)
    J = jacobian(pars,P,D,d)
    covmat = s2 * inv(transpose(J)*J)
    E = covmat[1:2,1:2]
    return x0, y0, E
end
export internochron
