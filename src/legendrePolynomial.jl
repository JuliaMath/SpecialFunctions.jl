function legendrePolynomial(n)
    #Storing the coefficients of the Legendre Polynomial up
    # to n degree in a n+1 x n+1 matrix. We need to store all the coefficients
    # because the loops are expression of the following recursive relation (Bonnet):
    #(n+1)*P_{n+1}(x) = (2n+1)*x*P_{n}(x) - n*P_{n-1}(x)
    c = zeros((n+1, n+1))
    c[1,1] = 1.0

    if n <= 0
        return c
    end
    c[2,2] = 1.0

    for i = 2:n
        for j=0:i-1
            c[i+1,j+1] = (- i + 1) * c[i-1,j+1] / i
        end
        for j=1:i
            c[i+1,j+1] = c[i+1,j+1] + (i+i - 1) * c[i,j] / i
        end
    end
    return c
end
