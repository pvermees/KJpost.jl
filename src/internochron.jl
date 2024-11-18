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
              xlim = [0,1.05*maximum([P./D;x0+sqrt(E[1,1])])],
              ylim = [0,1.05*maximum([d./D;y0+sqrt(E[2,2])])],
              legend = false,
              nsigma = 2,
              xlab = "P/D",
              ylab = "d/D",
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
    Plots.plot!([0,x0],[y0,0];seriescolor=:black,legend=legend)
    Plots.plot!(P./D,d./D;seriestype=:scatter,legend=legend,plot_options...)
    Plots.xlabel!(xlab)
    Plots.ylabel!(ylab)
    return p
end
function plot(x0,y0,E,
              P::AbstractVector,
              D::AbstractVector,
              d::AbstractVector,
              method::AbstractString;
              xlim = [0,1.05*maximum([P./D;x0+sqrt(E[1,1])])],
              ylim = [0,1.05*maximum([d./D;y0+sqrt(E[2,2])])],
              legend = false,
              nsigma = 2,
              plot_options...)
    if method=="U-Pb"
        t,st = x0y02t(x0,y0,E)
    else
        t,st = x02t(x0,sqrt(E[1,1]),method)
    end
    Pname,Dname,dname = Plasmatrace.getPDd(method)
    xlab = Pname * "/" * Dname
    ylab = dname * "/" * Dname
    p = PTpost.plot(x0,y0,E,P,D,d;
                    xlim=xlim,ylim=ylim,legend=legend,nsigma=nsigma,
                    xlab=xlab,ylab=ylab,plot_options...)
    sdig = 2
    tdig = ceil(Int,log10(t/st)) + sdig
    tstring = "t = " *
        string(round(t,sigdigits=tdig)) * "+/-" *
        string(round(st,sigdigits=sdig)) * "Ma"
    Plots.annotate!(Plots.xlims(p)[2],
                    Plots.ylims(p)[2],
                    Plots.text(tstring,:right))
    if method=="U-Pb"
        add_concordia_line(x0,y0)
    end
end
export plot

function add_concordia_line(x0,y0)
    L5, L8, U58 = UPb_helper()
    tmin = @. log(1+1/x0)/L8
    tmax = 4000
    t = range(tmin[1],tmax;step=50)
    x = @. 1/(exp(L8*t)-1)
    y = @. U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    Plots.plot!(x,y)
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

function UPb_helper()
    L5 = Plasmatrace._PT["lambda"]["U235-Pb207"][1]
    L8 = Plasmatrace._PT["lambda"]["U238-Pb206"][1]
    fUPb = Plasmatrace._PT["iratio"]["U-Pb"]
    U58 = fUPb.U235/fUPb.U238
    return L5, L8, U58
end

function x0y02t(x0::AbstractFloat,
                y0::AbstractFloat,
                E::Matrix)
    function misfit(t)
        L5, L8, U58 = UPb_helper()
        x = @. 1/(exp(L8*t)-1)
        y = @. U58*(exp(L5*t)-1)/(exp(L8*t)-1)
        yl = @. y0*(1-x/x0)
        return sum(@. (yl - y)^2)
    end
    fit = Optim.optimize(misfit,[1.0])
    t = Optim.minimizer(fit)[1]
    st = 0.1
    return t,st
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
export internochron
