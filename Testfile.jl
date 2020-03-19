using Pkg, PkgTemplates

t = Template(;
           user="ManuelStapper",
           license="MIT",
           authors=["Manuel Stapper <manuel.stapper@googlemail.com>"],
           dir="~/code",
           julia_version=v"0.7",
           ssh=false,
           plugins=[
               TravisCI(),
               Codecov(),
               Coveralls(),
               AppVeyor(),
           ],
       )

generate("Filtering", t)

using Distributions, Optim, LinearAlgebra, Plots, Random

Random.seed!(123)
x, y = kalman_sim(100, [1, 1], 1, [0.99 0; 0 0.6], [1 0; 0 1])
sm1 = kalman_smooth(y, [1, 1], 1, [0.99 0; 0 0.6], [1 0; 0 1])

plot(x[1, :], label = "True")
plot!(sm1[1][1, :], label = "Filter")
plot!(sm1[3][1, :], label = "Smooth")

kalman_EM(y, [1, 1], missing, missing, missing, 0.01)
