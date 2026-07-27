module deSolveDiffEqRCallExt

using deSolveDiffEq
using RCall

function __init__()
    deSolveDiffEq.rcall_adapter[] = RCall
    return nothing
end

end
