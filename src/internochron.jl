function get_pS(x0,y0,P,D,d)
    p = @. ((d*x0^2+d)*y0+D*x0^2-P*x0)/((P*x0+D)*y0^2+d*x0^2*y0+D*x0^2)
    S = @. (((d+P)*x0^2+(P+D)*x0+d+D)*y0^2+
            ((d+D)*x0^2+((-d)-P)*x0)*y0+(P+D)*x0^2)/((x0^2+1)*y0^2+x0^2)
    return p, S
end

function get_fitted_PDd(x0,y0,P,D,d)
    p, S = get_pS(x0,y0,P,D,d)
    x = @. x0*(1-p)
    y = @. y0*p
    z = @. 1+x+y
    Pf = @. S*x/z
    Df = @. S/z
    df = @. S*y/z
    return Pf, Df, df
end
export get_fitted_PDd

function jacobian(pars,P,D,d)
    ns = length(P)
    x0 = pars[1]
    y0 = pars[2]
    p, S = get_pS(x0,y0,P,D,d)
    x = @. x0*(1-p)
    y = @. y0*p
    z = @. 1+x+y
    J = zeros(2+2*ns,3*ns)
    J[1,1:ns] = @. S*(1-p)/z - (S*x0*(1-p)^2)/z^2                    # dPdx0
    J[1,ns+1:2*ns] = @. -S*(1-p)/z^2                                 # dDdx0
    J[1,2*ns+1:3*ns] = @. -S*p*(1-p)*y0/z^2                          # dddx0
    J[2,1:ns] = @. -S*p*(1-p)*x0/z^2                                 # dPdy0
    J[2,ns+1:2*ns] = @. -S*p/z^2                                     # dDdy0
    J[2,2*ns+1:3*ns] = @. S*p/z - (S*y0*p^2)/z^2                     # dddy0
    J[3:ns+2,1:ns] = diagm(@. x0*(1-p)/z)                            # dPdS
    J[3:ns+2,ns+1:2*ns] = diagm(@. 1/z)                              # dDdS
    J[3:ns+2,2*ns+1:3*ns] = diagm(@. p*y0/z)                         # dddS
    J[ns+3:end,1:ns] = diagm(@. -S*x0/z - S*(1-p)*x0*(y0-x0)/z^2)    # dPdp
    J[ns+3:end,ns+1:2*ns] = diagm(@. -S*(y0-x0)/z^2)                 # dDdp
    J[ns+3:end,2*ns+1:3*ns] = diagm(@. S*y0/z - S*p*y0*(y0-x0)/z^2) # dddp
    return J
end

function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector;
                      numerical::Bool=false)
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
        covmat = s2 * inv(J*transpose(J))
        E = covmat[1:2,1:2]
    end
    return x0, y0, E
end
function internochron(run::Vector{KJ.Sample},
                      method::AbstractString,
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple;
                      numerical::Bool=true)
    ns = length(run)
    out = DataFrame(t=fill(0.0,ns),st=fill(0.0,ns),
                    y0=fill(0.0,ns),sy0=fill(0.0,ns))
    for (i, samp) in enumerate(run)
        P, D, d = KJ.atomic(samp,channels,blank,pars)
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
    lambda = KJ._KJ["lambda"][method][1]
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
    L5 = KJ._KJ["lambda"]["U235-Pb207"][1]
    L8 = KJ._KJ["lambda"]["U238-Pb206"][1]
    fUPb = KJ._KJ["iratio"]["U-Pb"]
    U58 = fUPb.U235/fUPb.U238
    return L5, L8, U58
end
