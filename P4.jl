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

function gprime!(storage, c1)
    storage[1] = -1*(η/(η-1))*
    ((κ*(log(exp(c1[1])))^((η-1)/η)+(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((η-1)/η))^((η/(η-1))-1))*
    (((η-1)/η)*κ*(log(exp(c1[1])))^(((η-1)/η)-1)-((η-1)/η)*(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(((η-1)/η)-1)*(p1/p2))
end

function g2prime!(storage, c1)
    storage[1] = -1*(η/(η-1))*((η/(η-1))-1)*((κ*log(exp(c1[1]))^((η-1)/η))+(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((η-1)/η))^((η/(η-1))-2)*
    (((η-1)/η)*κ*log(exp(c1[1]))^(((η-1)/η)-1)-((η-1)/η)*(1-κ)*(p1/p2)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(((η-1)/η)-1))*
    (((η-1)/η)*κ*log((exp(c1[1])))^(((η-1)/η)-1)-((η-1)/η)*(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(((η-1)/η)-1)*(p1/p2))-
    (((η-1)/η)*(((η-1)/η)-1)*κ*log(exp(c1[1]))^(((η-1)/η)-2)+((η-1)/η)*(1-κ)*(((η-1)/η)-1)*log(exp((w-p1*log(exp(c1[1]))/p2)))^(((η-1)/η)-2)*(p1/p2)^2)*
    ((κ*log(exp(c1[1]))^((η-1)/η)+(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((η-1)/η))^((η/(η-1))-1))

end

function utility_maximizer(κ, η, p1, p2, w, initial_guess, method)
    if method == "Conjugate-Gradient"
        f(c1) = -1*(κ*c1[1]^((η-1)/η)+(1-κ)*((w-p1*c1[1])/p2)^((η-1)/η))^(η/(η-1))
        lower = [0.000001]
        upper = [w/p1]
        result = Optim.minimizer(optimize(f, fprime!, lower, upper, initial_guess, Fminbox(ConjugateGradient())))[1]
        return result
    elseif method == "Newton"
        g(c1) = -1*(κ*log(exp(c1[1]))^((η-1)/η)+(1-κ)*log(exp((w-p1*log(exp(c1[1]))/p2)))^((η-1)/η))^(η/(η-1))
        result = Optim.minimizer(optimize(g, gprime!, g2prime!, initial_guess, Newton()))[1]
        return result
    end
end

initial_guess = [2.0]
mtd = "Newton"
c1_star = utility_maximizer(κ, η, p1, p2, w, initial_guess, mtd)
c2_star = (w-p1*c1_star)/p2
