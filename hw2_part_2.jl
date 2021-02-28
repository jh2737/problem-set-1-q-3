# Step 0 set seed
# Step 1: Generate 20,000 * 20,0000 array of random numbers x
import Pkg; Pkg.add("Random");
using Random
Random.seed!(1234)
x = rand(Float64, (20000,20000))
num_rows = size(x)[1]
num_cols = size(x)[2]

# Step 2: Construct exp_cols which expoentiate elements of x column by column
function exp_cols(x)
    new_x = Array{Float64}(undef,num_rows,num_cols)
    for nth_col in 1:num_cols
        new_x[:,nth_col] = exp.(x[:,nth_col])
    end
    return new_x
end

# Step 3: Construct exp_rows which expoentiate elements of x row by row
function exp_rows(x)
    new_x = Array{Float64}(undef,num_rows,num_cols)
    for nth_row in 1:num_rows
        new_x[nth_row,:] = exp.(x[nth_row,:])
    end
    return new_x
end

# Step 4 time both functions
@time exp_cols(x)
@time exp_rows(x)
# exp_rows runs faster than exp_cols

# Step 5 compare with exponentiating the full Array
function exp_full(x)
    new_x = exp.(x)
    return new_x
end
@time exp_full(x) # exponentiating full array is the fastest
