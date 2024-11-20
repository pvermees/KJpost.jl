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

function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector;
                      numerical::Bool=true)
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
    if numerical
        J = ForwardDiff.jacobian(residuals,pars)
        E = s2 * inv(transpose(J) * J)
    else
        J = jacobian(pars,P,D,d)
        covmat = s2 * inv(transpose(J)*J)
        E = covmat[1:2,1:2]
    end
    return x0, y0, E
end
function internochron(run::Vector{Plasmatrace.Sample},
                      method::AbstractString,
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple;
                      numerical::Bool=true)
    ns = length(run)
    out = DataFrame(t=fill(0.0,ns),st=fill(0.0,ns),
                    y0=fill(0.0,ns),sy0=fill(0.0,ns))
    for (i, samp) in enumerate(run)
        P, D, d = Plasmatrace.atomic(samp,channels,blank,pars)
        x0, y0, E = internochron(P,D,d;numerical=numerical)
        out.t[i], out.st[i] = x0y02t(x0,y0,E,method)
        out.y0[i] = y0
        out.sy0[i] = sqrt(E[2,2])
    end
    return out
end
export internochron

function x0y02t(x0::AbstractFloat,
                y0::AbstractFloat,
                E::Matrix,
                method::AbstractString)
    if method == "U-Pb"
        return x0y02t(x0,y0,E)
    else
        return x02t(x0,sqrt(E[1,1]),method)
    end
end
function x02t(x0::AbstractFloat,
              sx0::AbstractFloat,
              method::AbstractString)
    lambda = Plasmatrace._PT["lambda"][method][1]
    R = 1/x0
    sR = R * sx0/x0
    t = log(1+R)/lambda
    st = (sR/(1+R))/lambda
    return t,st
end
function x0y02t(x0::AbstractFloat,
                y0::AbstractFloat,
                E::Matrix)
    L5, L8, U58 = UPb_helper()
    function misfit(par)
        t = par[1]
        x = 1/(exp(L8*t)-1)
        y = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
        yl = y0*(1-x/x0)
        return (yl - y)^2
    end
    init = log(1+1/x0)/L8
    fit = Optim.optimize(misfit,[init])
    t = Optim.minimizer(fit)[1]
    J = TWjacobian(t,x0,y0,E)
    covmat = J * E * transpose(J)
    st = sqrt(covmat[1,1])
    return t,st
end

function TWjacobian(t,x0,y0,E)
    L5, L8, U58 = UPb_helper()
    dfdx0 = -y0/x0^2
    dfdy0 = 1/x0 - exp(L8*t) + 1 
    dfdt = U58*L5*exp(L5*t) - L8*y0*exp(L8*t)
    dtdx0 = -dfdx0/dfdt
    dtdy0 = -dfdy0/dfdt
    J = hcat(dtdx0,dtdy0)
    return J
end

function UPb_helper()
    L5 = Plasmatrace._PT["lambda"]["U235-Pb207"][1]
    L8 = Plasmatrace._PT["lambda"]["U238-Pb206"][1]
    fUPb = Plasmatrace._PT["iratio"]["U-Pb"]
    U58 = fUPb.U235/fUPb.U238
    return L5, L8, U58
end
