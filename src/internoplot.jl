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
    t,st = x0y02t(x0,y0,E,method)
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
        xmax = Plots.xlims(p)[2]
        ymax = Plots.ylims(p)[2]
        add_concordia_line(xmax,ymax)
    end
end
export plot

function add_concordia_line(xmax,ymax)
    L5, L8, U58 = UPb_helper()
    tmin = log(1+1/xmax)/L8
    function Pb76misfit(par)
        t = par[1]
        ypred = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
        return (ymax - ypred)^2
    end
    Pb76fit = Optim.optimize(Pb76misfit,[4000.0])
    tmax = Optim.minimizer(Pb76fit)[1]
    t = range(tmin[1],tmax;step=50)
    x = @. 1/(exp(L8*t)-1)
    y = @. U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    Plots.plot!(x,y,linewidth=1.5,linecolor=:black)
end
