using Pkg
Pkg.add("Optim")
using Optim

# Parameters
κ = 0.5
η = 2.0
p1 = 1.0
p2 = 1.0
w = 10.0

function fprime!(storage, c1)
    storage[1] = -1*(η/(η-1))*
    ((κ*c1[1]^((η-1)/η)+(1-κ)*((w-p1*c1[1])/p2)^((η-1)/η))^((η/(η-1))-1))*
    (((η-1)/η)*κ*(c1[1])^(((η-1)/η)-1)-((η-1)/η)*(1-κ)*((w-p1*c1[1])/p2)^(((η-1)/η)-1)*(p1/p2))
end

function utility_maximizer(κ, η, p1, p2, w, initial_guess, method)
    if method == "Conjugate-Gradient"
        lower = [0.000001]
        upper = [w/p1]
        f(c1) = -1*(κ*c1[1]^((η-1)/η)+(1-κ)*((w-p1*c1[1])/p2)^((η-1)/η))^(η/(η-1))
        result = Optim.minimizer(optimize(f, fprime!, lower, upper, initial_guess, Fminbox(ConjugateGradient())))
        return result
    end
end

initial_guess = [1.0]
mtd = "Conjugate-Gradient"
outcome = utility_maximizer(κ, η, p1, p2, w, initial_guess, mtd)
