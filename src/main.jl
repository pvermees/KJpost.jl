function extend!(_PT::AbstractDict)

    if haskey(_PT["tree"]["extra"].action,"i")
        return # already added
    end
        
    extra_message = "i: Internal isochron\n" * _PT["tree"]["extra"].message
    extra_action = _PT["tree"]["extra"].action
    extra_action["i"] = "internochron"
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
        "p" => TUInternochronViewer!,
        "e" => TUInternochronExport
    )
    updateTree!(_PT["tree"],"internochron",
                message = internochron_message,
                help = internochron_help,
                action = internochron_action)

    view_message =
        "n: Next\n" *
        "p: Previous\n" *
        "g: Go to\n" *
        "t: Tabulate all the samples in the session\n" *
        "x: Exit\n" *
        "?: Help"
    view_help =
        "It is useful to view the outcome of the internal isochron " *
        "regression to ensure that the results are sensible and that the " *
        "algorithm did not get stuck in a local minimum."
    view_action = Dict(
        "n" => TUInternochron_next!,
        "p" => TUInternochron_previous!,
        "g" => "internogoto",
        "t" => Plasmatrace.TUItabulate
    )
    updateTree!(_PT["tree"],"internoview",
                message = view_message,
                help = view_help,
                action = view_action)

    goto_message =
        "Enter the number of the sample to plot " *
        "(? for help, x to exit):"
    goto_help = "Jump to a specific analysis."
    goto_action = TUInternochron_goto!
    updateTree!(_PT["tree"],"internogoto",
                message = goto_message,
                help = goto_help,
                action = goto_action)

    csv_message =
        "Enter the path and name of the .csv " *
        "file (? for help, x to exit):"
    csv_help =
        "Provide the file name with or without " *
        "the .csv extension."
    csv_action = internochron2csv
    updateTree!(_PT["tree"],"internochron2csv",
                message = csv_message,
                help = csv_help,
                action = csv_action)
    
end
export extend!

function TUInternochronViewer!(ctrl::AbstractDict)
    TUInternochron(ctrl)
    return "internoview"
end

function TUInternochron_next!(ctrl::AbstractDict)
    ctrl["i"] += 1
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    return TUInternochron(ctrl)
end

function TUInternochron_previous!(ctrl::AbstractDict)
    ctrl["i"] -= 1
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    return TUInternochron(ctrl)
end

function TUInternochron_goto!(ctrl::AbstractDict,
                              response::AbstractString)
    ctrl["i"] = parse(Int,response)
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    TUInternochron(ctrl)
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
    return nothing
end
export TUIinternochron

function TUInternochronExport(ctrl::AbstractDict)
    if "gui" in keys(ctrl)
        return GUIinternochronExport(ctrl)
    else
        return "internochron2csv"
    end
end
export TUInternochronExport

function GUIinternochronExport(ctrl::AbstractDict)
    save_dialog("Choose a csv file name",ctrl["gui"]) do fname
        @async internochron2csv(ctrl,fname)
    end
    return "xx"
end

function internochron2csv(ctrl::AbstractDict,
                          fname::AbstractString)
    print("I will export csv files")
end
