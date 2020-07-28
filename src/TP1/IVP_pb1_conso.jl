#
#
# Date:     juillet 2020
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# Intégration de l'equation differentiel du problème 1 avec minimisation de
# la consommation.
#
# ---------------

using Plots
using LinearAlgebra
using DifferentialEquations



include("lib/ode_rk4.jl")
include("lib/ode_euler.jl")
include("lib/ode_heun.jl")
include("lib/ode_runge.jl")

################################################################################
#############               Fonctions utiles                       #############
################################################################################

# -----------------------------------------------------------------------
# deuxieme membre de l'edo
# ------------------------
function hvfun_pb1_conso(t,z)   

    z2point = z[2]
    z1point = -z[1]
    if abs(z[2]) >= 1
        z1point=z1point+sign(z[2])
    end
    

    return [z1point z2point]'
end
#-----------------------------------------------------------------------
# meme fct mais ecrite differemment pour le solveur 
#-----------------------------------------------------------------------
function hvfun_ode_conso!(zpoint, z, p, t)
    z2point = z[2]
    z1point = -z[1]
    if abs(z[2]) >= 1
        z1point=z1point+sign(z[2])
    end
    zpoint[1] = z1point
    zpoint[2] = z2point

end

#
# Solution exacte de l'edo
# ------------------------
function exphvfun_pb1_conso(t,z0)

    if t>=1
        u = sign(z0[2])
    else
        u=0
    end
    zf = [z0[1]*exp(-t) + u*(1-exp(1-t))  z0[2]*exp(t)]

    return zf
end

#
# Fonction auxiliaires pour tracer les solutions
# --------------------
function plot_sol(plt, T, Y, c, label)

        plot!(T,Y[:,1], color=c, xlabel="t", ylabel="y_1(t)",label=label, subplot=1)
        plot!(T,Y[:,2], color=c, xlabel="t", ylabel="y_2(t)",label=label, subplot=2)
        plot!(Y[:,1],Y[:,2], color=c, xlabel="y_1(t)", ylabel="y_2(t)",label=label, subplot=3)

        display(plt)
end
#
# Fonction pause
#--------------------
pause(text) = (print(stdout, text); read(stdin, 1); nothing)


################################################################################
#############                 Début                                #############
################################################################################

x0 = 0
p0 = exp(-1)
z0 = [x0;  p0]
t0 = 0
tf = 2.1

N0 = [collect(120:60:1080); collect(1200:600:10800)]  # Valeurs de N pour les courbes d'ordre
N = 10



pyplot()
plt = Plots.plot(layout=(3))

#Euler
T, Z = ode_euler(hvfun_pb1_conso, [t0 tf], z0, N)
plot_sol(plt, T, Z, "magenta", "euler")

#################

#runge
T, X = ode_runge(hvfun_pb1_conso,[t0 tf],z0,N)
plot_sol(plt,T, X, "red", "runge")



################

#heun
T, Xheun = ode_heun(hvfun_pb1_conso,[t0 tf],z0,N)
plot_sol(plt, T, Xheun, "green", "heun")


##################

#rk4
T, Xrk4 = ode_rk4(hvfun_pb1_conso,[t0 tf],z0,N)
plot_sol(plt, T, Xrk4, "blue", "rk4")


#ODEProblem

RelTol = 1e-10
AbsTol = 1e-10

prob = ODEProblem(hvfun_ode_conso!, z0, (t0 ,tf))
sol = solve(prob,reltol=RelTol , abstol=AbsTol)

# extraction 
T = sol.t
XODE = sol.u
# Adaptation de la forme de X
X = XODE[1]'
for i in 2:length(XODE)
   global X = [X;XODE[i]']
end


plot_sol(plt, T, X, "cyan", "ODEProblem")


################################################################################
######                         ORDRE                                     #######
################################################################################

pause("tapez entrée pour voir le graphique des ordres")
plt = Plots.plot(layout=(1,2))


#Courbes d'ordre
N=N0
zf = exphvfun_pb1_conso(tf,z0)

#erreur [ordre]
err1=zeros(length(N),2)
err2=err1
err3=err1
err4=err1
err5 = err4

#Euler
nfe= N
for i = 1:length(N)
    T, Z = ode_euler(hvfun_pb1_conso, [t0 tf], z0, N[i])
    err1[i,:] = log10.(abs.( Z[end,:]' - zf))
end


Plots.plot!(log10.(nfe),err1[:,1],color="magenta" ,xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="Euler")
Plots.plot!(log10.(nfe),err1[:,2],color="magenta" , xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="Euler")


#Runge
for i = 1:length(N)
    T, Z=ode_runge(hvfun_pb1_conso, [t0 tf], z0, N[i]/2)
    err2[i,:] = log10.(abs.( Z[end,:]' - zf))
end


Plots.plot!(log10.(nfe),err2[:,1],color="red" ,xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="runge")
Plots.plot!(log10.(nfe),err2[:,2],color="red" , xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="runge")


#Heun
for i = 1:length(N)
    T,Z=ode_heun(hvfun_pb1_conso, [t0 tf], z0, N[i]/3)
    err3[i,:] = log10.(abs.( Z[end,:]' - zf))
end


Plots.plot!(log10.(nfe),err3[:,1],color="green" ,xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="heun")
Plots.plot!(log10.(nfe),err3[:,2],color="green" , xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="heun")



#RK4
for i = 1:length(N)
    T, Z=ode_rk4(hvfun_pb1_conso, [t0 tf], z0, N[i]/4)
    err4[i,:] = log10.(abs.( Z[end,:]' - zf))
end



Plots.plot!(log10.(nfe),err4[:,1],color="blue" ,xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="rk4")
Plots.plot!(log10.(nfe),err4[:,2],color="blue" , xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="rk4")

#ODEProblem
for i = 1:length(N)

    prob = ODEProblem(hvfun_ode_conso!, z0, (t0 ,tf))
    sol = solve(prob,reltol=RelTol , abstol=AbsTol)

    # extraction 
    T = sol.t
    XODE = sol.u
    # Adaptation de la forme de X
    err5[i,:] = log10.(abs.(XODE[end]' - zf))
end

Plots.plot!(log10.(nfe),err5[:,1],color="cyan" ,xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_1)",subplot=1,label="ODEProblem")
Plots.plot!(log10.(nfe),err5[:,2],color="cyan" , xlabel="log_{10}(fe)", ylabel="log_{10}(erreur pour y_2)",subplot=2,label="ODEProblem")





