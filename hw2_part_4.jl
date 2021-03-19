
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

function g_first!(storage, c1)
    storage[1] = -1*(e2)*
    ((kappa*(log(exp(c1[1])))^(e1)+(1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(e1))^((e2)-1))*
    ((e1)*kappa*(log(exp(c1[1])))^((e1)-1)-(e1)*(1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((e1)-1)*(p1/p2))
end

function g_second!(storage, c1)
    storage[1] = -1*(e2)*((e2)-1)*((kappa*log(exp(c1[1]))^(e1))+(1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(e1))^((e2)-2)*
    ((e1)*kappa*log(exp(c1[1]))^((e1)-1)-(e1)*(1-kappa)*(p1/p2)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((e1)-1))*
    ((e1)*kappa*log((exp(c1[1])))^((e1)-1)-(e1)*(1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((e1)-1)*(p1/p2))-
    ((e1)*((e1)-1)*kappa*log(exp(c1[1]))^((e1)-2)+(e1)*(1-kappa)*((e1)-1)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((e1)-2)*(p1/p2)^2)*
    ((kappa*log(exp(c1[1]))^(e1)+(1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(e1))^((e2)-1))

end


function utility_maximizer(kappa, eta, p1, p2, w, initial_guess, solver)
    if solver == "Newton"
        g(c1) =  -(kappa*log(exp(c1[1])))^e1 - ((1-kappa)*log(exp((w-p1*log(exp(c1[1]))/p2)))^e1)^e2
        output = Optim.minimizer(optimize(g, g_first!, g_second!, initial_guess, Newton()))
    elseif solver == "Conjugate gradient"
        lower = [1e-3]
        upper = [w/p1]
        f(c1) = -1*(kappa*c1[1]^e1+(1-kappa)*((w-p1*c1[1])/p2)^e1)^e2
        output = Optim.minimizer(optimize(f, f_first!, lower, upper, initial_guess, Fminbox(ConjugateGradient())))
    end
    return output
end

utility_maximizer(kappa, eta, p1, p2, w, initial_guess, "Newton")
