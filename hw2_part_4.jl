
# Utility function without constraint
using Optim

kappa = 0.5
eta = 2.0
p1 = 1.0
p2 = 1.0
w = 10.0
e1 = (eta-1)/eta
e2 = eta/(eta-1)
initial_guess = [1.0]

function f_first!(storage, c1)
    storage[1] = -1*e2*
    ((kappa*c1[1]^(e1)+(1-kappa)*((w-p1*c1[1])/p2)^(e1))^(e2-1))*
    (e1*kappa*(c1[1])^((e1)-1)-e1*(1-kappa)*((w-p1*c1[1])/p2)^(e1-1)*(p1/p2))
end


function utility_maximizer(kappa, eta, p1, p2, w, initial_guess, solver)

    f(c1) = -1*(kappa*c1[1]^e1+(1-kappa)*((w-p1*c1[1])/p2)^e1)^e2
    lower = [1e-3]
    upper = [w/p1]

    if solver == "Conjugate gradient"
        output = Optim.minimizer(optimize(f, f_first!, lower, upper, initial_guess, Fminbox(ConjugateGradient())))
    end
    return output
end

utility_maximizer(kappa, eta, p1, p2, w, initial_guess, "Conjugate gradient")
