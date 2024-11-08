function misfit(par::AbstractVector,
                P::AbstractVector,
                D::AbstractVector,
                d::AbstractVector)
    x0 = par[1]
    y0 = par[2]
    p = @. ((d*x0^2+d)*y0+D*x0^2-P*x0)/((P*x0+D)*y0^2+d*x0^2*y0+D*x0^2)
    Df = @. ((P*x0+D)*y0^2+d*x0^2*y0+D*x0^2)/((x0^2+1)*y0^2+x0^2)
    SS = @. (d-Df*p*y0)^2+(P-Df*(1-p)*x0)^2+(D-Df)^2
    return sum(S)
end

function internochron(P::AbstractVector,
                      D::AbstractVector,
                      d::AbstractVector)
    X = P./D
    Y = d./D
    y0i, slope = linear_regression(X,Y)
    x0i = -y0i/slope
    fit = Optim.optimize(misfit,[x0,y0])
    [x0,y0] = Optim.minimizer(fit)
    return x0, y0
end

function linear_regression(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)
    numerator = sum((x .- x_mean) .* (y .- y_mean))
    denominator = sum((x .- x_mean).^2)
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    return intercept, slope
end
