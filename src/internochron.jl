function get_fitted_PDd(x0,y0,P,D,d)
    p = @. ((d*x0^2+d)*y0+D*x0^2-P*x0)/((P*x0+D)*y0^2+d*y0*x0^2+D*x0^2)
    Df = @. ((P*x0+D)*y0^2+d*y0*x0^2+D*x0^2)/((x0^2+1)*y0^2+x0^2)
    Pf = @. Df*x0*(1-p)
    df = @. Df*y0*p
    return Pf, Df, df
end
export get_fitted_PDd

function plot(x0::AbstractFloat,
              y0::AbstractFloat,
              P::AbstractVector,
              D::AbstractVector,
              d::AbstractVector)
    p = Plots.plot(P./D,d./D,seriestype=:scatter,legend=false)
    Plots.plot!([0,x0],[y0,0],legend=false)
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
        return [P-Pf;D-Df;d-df]
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
    J = ForwardDiff.jacobian(residuals,pars)
    E = s2 * inv(transpose(J) * J)
    sx0 = sqrt(E[1,1])
    sy0 = sqrt(E[2,2])
    rho = E[1,2]/sqrt(E[1,1]*E[2,2])
    return x0, sx0, y0, sy0, rho
end
export internochron
