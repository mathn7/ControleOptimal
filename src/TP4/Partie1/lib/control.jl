function control(t, z, par)
    #-------------------------------------------------------------------------------------------
    #
    #    control()
    #
    #    Description
    #
    #        Computes the control.
    #
    #-------------------------------------------------------------------------------------------
    #
    #    Usage
    #
    #        u = control(t, z, par)
    #
    #    Inputs
    #
    #        t    -  real        : time
    #        z    -  real vector : state & costate
    #        par  -  real vector : parameters, par=[] if no parameters
    #
    #    Outputs
    #
    #        u   -  real vector : control()
    #
    #-------------------------------------------------------------------------------------------
    
    # par = [t0; tf; x0; xf, epsilon]'
    #
    n  = length(z)/2
    x  = z[1:n]
    p  = z[n+1:2*n]
    
    epsilon = par[5]
    w = -1 + abs(p)
    if p!=0
        u = (-2*epsilon*sign(p))/(w-2*epsilon-sqrt(w^2+4*epsilon^2))
    else
        u = (-2*epsilon)/(-1-2*epsilon-sqrt(1+4*epsilon^2))
        
    end

    return u
end
