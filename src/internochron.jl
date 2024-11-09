function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector)
    function misfit(par::AbstractVector)
        x0 = par[1]
        y0 = par[2]
        p = @. ((d*x0^2+d)*y0+D*x0^2-P*x0)/((P*x0+D)*y0^2+d*x0^2*y0+D*x0^2)
        Df = @. ((P*x0+D)*y0^2+d*x0^2*y0+D*x0^2)/((x0^2+1)*y0^2+x0^2)
        SS = @. (d-Df*p*y0)^2+(P-Df*(1-p)*x0)^2+(D-Df)^2
        return sum(SS)
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
    return x0, y0
end
export internochron
