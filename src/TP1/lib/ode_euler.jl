@doc doc"""

# Syntaxe
```julia
T, X = ode_euler(f,t0tf,y0,N)```
# Entrée :
   * **f**    : Function -second member of the ode whith the interface -xpoint = f(t, x) t    - real     : time, x = vector of R^n with the same dimension of x0
   * **t0tf** : intial and final time  [t0,tf]
   * **y0**   : Array{Float64,1}      : initial point
   * **N**    : Int number of steps (>1)

# Sortie:
   * **T**    : Array{Float64,1} real(N+1,1)  : vector of times
   * **X**    : Array{Float64,2}    Matrix of solution

"""


function  ode_euler(f::Function,t0tf,x0,N)

    N = Int(N)
    n = length(x0)
    T = zeros(N+1,1)
    X = zeros(N+1,n)

    pas = (t0tf[2]-t0tf[1])/N
    T[1] = 0
    X[1,:] = x0

    for i = 2:N+1
        X[i,:] = X[i-1,:] + pas*f(T[i-1],X[i-1,:])
        T[i] = (i-1)*pas
    end

    return T, X
end

