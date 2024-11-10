function extend!(_PT::AbstractDict)

    if !haskey(_PT["tree"]["extra"].action,"i")
        
        extra_message = "i: Internal isochron\n" * _PT["tree"]["extra"].message
        extra_action = _PT["tree"]["extra"].action
        extra_action["i"] = "internal_isochron"
        updateTree!(_PT["tree"],"extra",
                    message = extra_message,
                    action = extra_action)
        internochron_message =
            "p: Plot\n" *
            "e: Export\n" *
            "x: Exit\n" *
            "?: Help"
        internochron_help =
            "Plot the internal isochron for a single laser spot " *
            "or export all the spots to a .csv file."
        internochron_action = Dict(
            "p" => TUInternochron,
            "e" => TUInternochronExport
        )
        updateTree!(_PT["tree"],"internal_isochron",
                    message = internochron_message,
                    help = internochron_help,
                    action = internochron_action)
        
        view_message = 
            "n: Next\n"*
            "p: Previous\n"*
            "g: Go to\n"*
            "t: Tabulate all the samples in the session\n"*
            "x: Exit\n"*
            "?: Help"
        view_help =
            "It is useful to view the outcome of the internal isochron "*
            "regression to ensure that the results are sensible and that the "*
            "algorithm did not get stuck in a local minimum."
        view_action = Dict(
            "n" => TUInternochron_next!,
            "p" => TUInternochron_previous!,
            "g" => "internochron_goto",
            "t" => TUItabulate
        )
        
    end

end
export extend!

function TUInternochron_next!(ctrl::AbstractDict)
    return "x"
end

function TUInternochron_previous!(ctrl::AbstractDict)
    return "x"
end

function updateTree!(tree::AbstractDict,
                     key::AbstractString;
                     message=tree[key].message,
                     help=tree[key].help,
                     action=tree[key].action)
    tree[key] = (message=message,help=help,action=action)
end

function TUInternochron(ctrl::AbstractDict)
    P, D, d = Plasmatrace.atomic(ctrl["run"][ctrl["i"]],
                                 ctrl["channels"],
                                 ctrl["blank"],
                                 ctrl["par"])
    x0, y0, E = internochron(P,D,d)
    p = plot(x0,y0,E,P,D,d)
    display(p)
    return "x"
end
export TUIinternochron

function TUInternochronExport(ctrl::AbstractDict)
end
