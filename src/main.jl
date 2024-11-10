function extend!(_PT::AbstractDict)

    if !haskey(_PT["tree"]["extra"].action,"i")
        
        message1 = "i: Internal isochron\n" * _PT["tree"]["extra"].message
        action1 = _PT["tree"]["extra"].action
        action1["i"] = "internal_isochron"
        updateTree!(_PT["tree"],"extra",
                    message = message1,
                    action = action1)
        message2 =
            "p: Plot\n" *
            "e: Export\n" *
            "x: Exit\n" *
            "?: Help"
        help2 =
            "Plot the internal isochron for a single laser spot " *
            "or export all the spots to a .csv file."
        action2 = Dict(
            "p" => TUInternochron,
            "e" => TUInternochronExport
        )
        updateTree!(_PT["tree"],"internal_isochron",
                    message = message2,
                    help = help2,
                    action = action2)
        
    end

end
export extend!

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
