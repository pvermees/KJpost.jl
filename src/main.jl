using Infiltrator

function extend!(_PT::AbstractDict)

    message = "i. Internal isochron\n" * _PT["tree"]["extra"].message
    action = _PT["tree"]["extra"].action
    action["i"] = TUIintchron
    
    updateTree!(_PT["tree"],"extra",
                message = message,
                action = action)

end
export extend!

function updateTree!(tree::AbstractDict,
                     key::AbstractString;
                     message=tree[key].message,
                     help=tree[key].help,
                     action=tree[key].action)
    tree[key] = (message=message,help=help,action=action)
end

function TUIintchron(ctrl::AbstractDict)
    print("I will carry out internal isochron regression")
end
export TUIintchron
