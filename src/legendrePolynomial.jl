function legendrePolynomial(x,n)
    #Evaluating the legendre polynomial of degree n in the point x.
    sum = 0
    for k = 0:n
        sum += (factorial(n)/(factorial(n-k)*factorial(k)))^2*(x-1)^(n-k)*(x+1)^k
    end
    evaluation = 1/(2^n)*sum
    return evaluation
end
