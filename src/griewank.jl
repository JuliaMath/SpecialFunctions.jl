using Base.Math: @horner

# Computes Griewank Scores
# Griewank Function is used in testing of optimization. 
# Defined in: https://en.wikipedia.org/wiki/Griewank_function

"""
griewank(x::Float64)

Computes first order griewank score for particular value of x.
Please use griewankm(A::Array{Float64, 1}) for first oder griewank calculations.
"""
function griewank(x::Float64)
    g = 1 + 1/4000*(x^2) - cos(x)
end

"""
griewank(A::Array{Float64, 1})

Computes first order griewank function for the values in the supplied array.
Returns an array with first order griewank scores.
These scores may be easily plotted, for example plot(griewank(A)).
"""
function griewankm(A::Array{Float64, 1})
    s = size(A)[1]
    resultant = Float64[]
    for i in 1:s
        push!(resultant, griewank(A[i]))
    end
    resultant
end

"""
griewank(A::Array{Float64, 1}, n::Int64)

Computes griewank score of n-th order. 
Please use griewankmn(A::Array{Float64, 2}) to calculate griewank scores of n-th order.
"""
function griewankn(A, n::Int64)
    sum = 0
    prod = 1
    for i in 1:n
        sum = sum + A[i]^2
    end
    for i in 1:n
       prod = prod * cos(A[i]/sqrt(i))
    end
    r = 1 + 1/4000*sum - prod
end

"""
griewank(A::Array{Float64, 2}) 

Computes griewank scores for n-th order. 
Accepts a matrix of order MxN and returns a vector with n-th order griewank scores.
"""
function griewankmn(A::Array{Float64, 2})
    M = size(A)[1]
    N = size(A)[2]
    resultant = Float64[]
    for i in 1:M
        push!(resultant, griewankn(view(A,i,:), N))
    end
    resultant
end
