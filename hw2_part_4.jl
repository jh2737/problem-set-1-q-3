
# Utility function without constraint
using Optim, ForwardDiff
function utility_maximizer(kappa, eta, p1, p2, w, c1, solver)
    f(c1) = (kappa*c1^(1-1/eta) + (1-kappa)*((w-p1*c1)/p2)^(1-1/eta))^(eta/(eta-1))
    if solver == "Newton"
        output = optimize(f, ForwardDiff.gradient(f, c1), ForwardDiff.hessian(f, c1))
    elseif solver == "Conjugate gradient"
        output = optimize(f, ForwardDiff.gradient(f, c1), ConjugateGradient())
    end
    return output
end
###
