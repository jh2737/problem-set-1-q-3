
# Utility function without constraint
using Optim, ForwardDiff
function utility_maximizer(kappa, eta, p1, p2, w, c1, solver)
    e1 = (eta-1)/eta
    e2 = eta/(eta-1)
    f(c1) = kappa*c1^e1+(1-kappa)*((w-p1*c1)/p2)^e1
    u(c1) = -f(c1)^e2
    # first derivative function
    # u'(c1) = -e2 * f(c1)^(e2-1) * (e1*kappa*c1^(e1-1) + e1*(1-kappa)*((w-p1*c1)/p2)^(e1-1))

    function u!(storage, c1)
        storage[1] = -e2 * f(c1)^(e2-1) * (e1*kappa*c1^(e1-1) + e1*(1-kappa)*((w-p1*c1)/p2)^(e1-1))
    end

    # second derivative function
    # second derivative
    # set h(c1) = e1*kappa*c1^(e1-1) + e1*(1-kappa)*((w-p1*c1)/p2)^(e1-1)
    # u''(c1) = -e2*(e2-1)*f(c1)^(e2-2)*h(c1)*h(c1) - e2*f(c1)^(e2-1)*(e1*(e1-1)*kappa*c1^(e1-2)+e1*(e1-1)*(1-kappa)*(w-p1*c1)/p2)^(e1-2))
    h(c1) = e1*kappa*c1^(e1-1) + e1*(1-kappa)*((w-p1*c1)/p2)^(e1-1)
    function uu!(storage, c1)
        storage[1] = -e2*(e2-1)*f(c1)^(e2-2)*h(c1)*h(c1) -
        e2*f(c1)^(e2-1)*(e1*(e1-1)*kappa*c1^(e1-2)+e1*(e1-1)*(1-kappa)*(w-p1*c1)/p2)^(e1-2))
    end
    if solver == "Newton"
        output = optimize(u, u!, uu!, 0)
    elseif solver == "Conjugate gradient"
        output = optimize(c1->f(first(c1)),[initial_guess], ConjugateGradient())
    end
    return output
end
