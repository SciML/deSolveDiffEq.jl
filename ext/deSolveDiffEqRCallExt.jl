module deSolveDiffEqRCallExt

using deSolveDiffEq
using RCall

function __init__()
    deSolveDiffEq.r_adapter[] = RCall
    return nothing
end

end
