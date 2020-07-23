#
#
# Auteurs:  Joseph GERGAUD
# Date:     juillet 2020
# Adresse:  INP-ENSEEIHT-IRIT-UMR CNRS 5055
#           2; rue Camichel 31071 Toulouse FRANCE
# Email:    gergaud@enseeiht.fr
#***************************************************************************
#
# Intégration de l'equation differentiel du problème 1 avec minimisation de
# l'énergie.
#
# ---------------

using Plots
using LinearAlgebra
using PyPlot

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
function hvfun_pb1_energie(t,z)
    zpoint = [-z[1]+z[2]  z[2]]'
    return zpoint
end

#
# Solution exacte de l'edo
# ------------------------
function exphvfun_pb1_energie(t,z0)
    #zf = [(z0[2]/2)*exp(t)+(z0[1]-z0[2]/2)*exp(-t); z0[2]*exp(t)]
    zf = zeros(1,2)
    zf[1] = (z0[2]/2)*exp(t)+(z0[1]-z0[2]/2)*exp(-t)
    zf[2] = z0[2]*exp(t)
    return zf
end

#
# Fonction auxiliaires pour tracer les solutions
# --------------------
function plot_sol(plt, T, Y, c, label, Y1Y2)
    if Y1Y2
        plot!(T,Y[:,1], color=c, xlabel="t", ylabel="y_1(t)",label=label, subplot=1)
        plot!(T,Y[:,2], color=c, xlabel="t", ylabel="y_2(t)",label=label, subplot=2)
        plot!(Y[:,1],Y[:,2], color=c, xlabel="y_1(t)", ylabel="y_2(t)",label=label, subplot=3)
    else
       plot!(T,Y[:,1], color=c, xlabel="log10(fe)", ylabel="log10(erreur pour y1)",label=label, subplot=1)
       plot!(T,Y[:,2], color=c, xlabel="log10(fe)", ylabel="log10(erreur pour y2)",label=label, subplot=2)
    end
    display(plt)
end

################################################################################
#############                 Début                                #############
################################################################################

x0 = 0
p0 = exp(-1)
z0 = [x0;  p0]
t0 = 0
tf = 2.1
N0 = [120:60:1080 1200:600:10800] # Valeurs de N pour les courbes d'ordre
N = 10


closeall()
pyplot()
plt = Plots.plot(layout=(1,3))

T, Z = ode_euler(hvfun_pb1_energie, [t0 tf], z0, N)
plot_sol(plt, T, Z, "magenta", "euler", true)

#################


T, X = ode_runge(hvfun_pb1_energie,[t0 tf],z0,N)
plot_sol(plt,T, X, "red", "runge", true)



################


T, Xheun = ode_heun(hvfun_pb1_energie,[t0 tf],z0,N)
plot_sol(plt, T, Xheun, "green", "heun", true)


##################

T, Xrk4 = ode_rk4(hvfun_pb1_energie,[t0 tf],z0,N)
plot_sol(plt, T, Xrk4, "blue", "rk4", true)



################################################################################
######                         ORDRE                                     #######
################################################################################

#Courbes d'ordre
N=N0
zf = exphvfun_pb1_energie(tf,z0)

#erreur [ordre]
err1=zeros(length(N),2)
err2=err1
err3=err1
err4=err1
err5 = err4

#Euler
nfe= N
for i = 1:length(N)
    T, Z = ode_euler(hvfun_pb1_energie, [t0 tf], z0, N[i])
    err1[i,:] = log10.(abs.( Z[end,:]' - zf))
end
#=
#Runge
for i = 1:length(N)
    T, Z=ode_runge(hvfun_pb1_energie, [t0 tf], z0, N[i]/2)
    err2[i,:] = log10.(abs.( Z[end,:]' - zf))
end

#Heun
for i = 1:length(N)
    T,Z=ode_heun(hvfun_pb1_energie, [t0 tf], z0, N[i]/3)
    err3[i,:] = log10.(abs.( Z[end,:]' - zf))
end

# #RK4
for i = 1:length(N)
    T, Z=ode_rk4(hvfun_pb1_energie, [t0 tf], z0, N[i]/4)
    err4[i,:] = log10.(abs.( Z[end,:]' - zf))
end

pyplot()
plt = Plots.plot(layout=(1,2))

plot_sol(plt, log10.(nfe), err1, "magenta", "euler", false )
#=
plot_sol(plt, log10.(nfe), err2, "red", "runge", false )
plot_sol(plt, log10.(nfe), err3, "green", "heun", false )
plot_sol(plt, log10.(nfe), err4, "blue", "rk4", false )
=#
