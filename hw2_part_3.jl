
using ForwardDiff
function cournot_newton(eta, c1, c2, guess, tolerance) # note that guess is a 1*2 array
    profit_1(guess) = (sum(guess))^(-1/eta)*guess[1] - 0.5*c1*guess[1]^2
    profit_2(guess) = (sum(guess))^(-1/eta)*guess[2] - 0.5*c2*guess[2]^2

    diff = 1e3
    while diff > tolerance #here diff is the distance matrix
        println("Now the guess for q is $guess.")
        f1 = guess -> ForwardDiff.gradient(profit_1,guess) #calculate FOC w.r.t q1
        f2 = guess -> ForwardDiff.gradient(profit_2,guess) #calculate FOC w.r.t q2
        f1_prime = guess -> ForwardDiff.hessian(profit_1,guess) #calculate SOC w.r.t q1
        f2_prime = guess -> ForwardDiff.hessian(profit_2,guess) #calculate SOC w.r.t q2
        new_guess_q1 = guess[1] - f1(guess)[1] / f1_prime(guess)[1,1]
        new_guess_q2 = guess[2] - f2(guess)[2] / f2_prime(guess)[2,2]
        diff = sqrt((guess[1]-new_guess_q1)^2 + (guess[2]-new_guess_q2)^2)
        guess = [new_guess_q1 new_guess_q2]
    end
    println("The optimal quantity is $guess")
end
