include("control.jl")
function hvfun(t, z, par)
    #-------------------------------------------------------------------------------------------
    #
    #    hvfun()
    #
    #    Description
    #
    #        Computes the Hamiltonian vector field associated to H.
    #
    #-------------------------------------------------------------------------------------------
    #
    #    Usage
    #
    #        hv = hvfun(t, z, par)
    #
    #    Inputs
    #
    #        t    -  real        : time
    #        z    -  real vector : state & costate
    #        par  -  real vector : parameters, par=[] if no parameters
    #
    #    Outputs
    #
    #        hv   -  real vector : hamiltonian vector field at time t
    #
    #-------------------------------------------------------------------------------------------
    n  = length(z)/2
    x  = z[1:n]
    p  = z[n+1:2*n]
    u = control[t,z,par]
    
    hv = [-x+u  p]'

    return hv
    end
    