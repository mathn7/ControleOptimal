using DifferentialEquations

function exphvfun(tspan, z0, options, par)
    #-------------------------------------------------------------------------------------------
    #
    #    exphvfun()
    #
    #    Description
    #
    #        Computes the chronological exponential of the Hamiltonian vector field hv
    #        defined by h.
    #
    #-------------------------------------------------------------------------------------------
    #
    #    Matlab Usage
    #
    #        [tout, z, flag] = exphvfun(tspan, z0, options, par)
    #
    #    Inputs
    #
    #        tspan   - real row vector of dimension 2   : tspan = [t0 tf]
    #        z0      - real vector                      : initial flow()
    #        options - struct vector                    : odeset options
    #        par     - real vector                      : parameters, par=[] if no parameters
    #
    #    Outputs
    #
    #        tout    - real row vector                  : time at each integration step
    #        z       - real matrix                      : z[:,i] : flow at tout[i]
    #
    #-------------------------------------------------------------------------------------------
    

    prob = ODEProblem(hvfun, z0, (t0 ,tf))

    sol = solve(prb, reltol= options[1], abstol = options[2])
    tout = sol.t
    x = sol.u 
    # Adaptation de la forme de X
    X = x[1]'
    for i in 2:length(x)
        global z = [z;x[i]']
    end
    return tout, z

end

