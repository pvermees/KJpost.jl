function get_fitted_PDd(x0::AbstractFloat,
                        y0::AbstractFloat,
                        P::AbstractVector,
                        D::AbstractVector,
                        d::AbstractVector)
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

function hessian(x0,y0,Pm,Dm,dm)
    p = @. ((dm*x0^2+dm)*y0+Dm*x0^2-Pm*x0)/((Pm*x0+Dm)*y0^2+dm*x0^2*y0+Dm*x0^2)
    D = @. ((Pm*x0+Dm)*y0^2+dm*x0^2*y0+Dm*x0^2)/((x0^2+1)*y0^2+x0^2)
    d2SSdx02 = sum(@. 2*D^2*(1-p)^2)
    d2SSdx0dy0 = 0
    d2SSdx0dp = @. 2*D*(Pm-D*(1-p)*x0)-2*D^2*(1-p)*x0
    d2SSdx0dD = @. 2*D*(1-p)^2*x0-2*(1-p)*(Pm-D*(1-p)*x0)
    d2SSdy0dx0 = 0
    d2SSdy02 = sum(@. 2*D^2*p^2)
    d2SSdy0dp = @. 2*D^2*p*y0-2*D*(dm-D*p*y0)
    d2SSdy0dD = @. 2*D*p^2*y0-2*p*(dm-D*p*y0)
    d2SSdpdx0 = @. 2*D*(Pm-D*(1-p)*x0)-2*D^2*(1-p)*x0
    d2SSdpdy0 = @. 2*D^2*p*y0-2*D*(dm-D*p*y0)
    d2SSdp2 = @. 2*D^2*y0^2+2*D^2*x0^2
    d2SSdpdD = @. 2*D*p*y0^2-2*y0*(dm-D*p*y0)-2*D*(1-p)*x0^2+2*x0*(Pm-D*(1-p)*x0)
    d2SSdDdx0 = @. 2*D*(1-p)^2*x0-2*(1-p)*(Pm-D*(1-p)*x0)
    d2SSdDdy0 = @. 2*D*p^2*y0-2*p*(dm-D*p*y0)
    d2SSdDdp = @. 2*D*p*y0^2-2*y0*(dm-D*p*y0)-2*D*(1-p)*x0^2+2*x0*(Pm-D*(1-p)*x0)
    d2SSdD2 = @. 2*p^2*y0^2+2*(1-p)^2*x0^2+2
    ns = length(Pm)
    ix0 = 1
    iy0 = 2
    ip = 3:2+ns
    iD = 3+ns:2+2*ns
    H = zeros(2*(ns+1),2*(ns+1))
    H[ix0,ix0] = d2SSdx02
    H[ix0,iy0] = d2SSdx0dy0
    H[ix0,ip] = d2SSdx0dp
    H[ix0,iD] = d2SSdx0dD
    H[iy0,ix0] = d2SSdy0dx0
    H[iy0,iy0] = d2SSdy02
    H[iy0,ip] = d2SSdy0dp
    H[iy0,iD] = d2SSdy0dD
    H[ip,ix0] = d2SSdpdx0
    H[ip,iy0] = d2SSdpdy0
    H[ip,ip] = diagm(d2SSdp2)
    H[ip,iD] = diagm(d2SSdpdD)
    H[iD,ix0] = d2SSdDdx0
    H[iD,iy0] = d2SSdDdy0
    H[iD,ip] = diagm(d2SSdDdp)
    H[iD,iD] = diagm(d2SSdD2)
    return H
end

function jacobian(x0,y0,Pm,Dm,dm)
    p = @. ((dm*x0^2+dm)*y0+Dm*x0^2-Pm*x0)/((Pm*x0+Dm)*y0^2+dm*x0^2*y0+Dm*x0^2)
    D = @. ((Pm*x0+Dm)*y0^2+dm*x0^2*y0+Dm*x0^2)/((x0^2+1)*y0^2+x0^2)
    dSSdx0 = sum(@. -2*D*(1-p)*(Pm-D*(1-p)*x0))
    dSSdy0 = sum(@. -2*D*p*(dm-D*p*y0))
    dSSdp = 2*D*x0*(Pm-D*(1-p)*x0)-2*D*y0*(dm-D*p*y0)
    dSSdD = (-2*p*y0*(dm-D*p*y0))-2*(1-p)*x0*(Pm-D*(1-p)*x0)-2*(Dm-D)
    ns = length(Pm)
    return [dSSdx0;dSSdy0;dSSdp;dSSdD]
end

function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector)
    function misfit(par::AbstractVector)
        Pf, Df, df = get_fitted_PDd(par[1],par[2],P,D,d)
        SS = @. (df-d)^2+(Pf-P)^2+(Df-D)^2
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
    J = jacobian(x0,y0,P,D,d)
    H = hessian(x0,y0,P,D,d)
    E = inv(H)
    return x0, y0
end
export internochron
