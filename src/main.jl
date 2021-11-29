#=
main:
- Julia version: 
- Author: sutymate
- Date: 2021-11-27
=#

include("../src/inputs.jl")

function apx(index, row, parameters, b)
    # solves for x_index on row_index with parameters {x_1, ..., x_MATRIX_DIMENSION}
    # lhs is x_index
    # rhs are all variables {x_1, ..., x_MATRIX_DIMENSION} moved from equation Ax = b, and b_index
    sum = Float64(0.0)
    c = row[index]
    # x_index has constant parameter 0, because it is on the left side of equation
    row[index] = Float64(0.0)
    for i in 1:MATRIX_DIMENSION
        sum += -1*row[i]*parameters[i]
    end
    sum += b[index]
    # prevent zero division
    if c != 0.0 sum /= c end
    return sum
end


function jacobi(R, D_, parameters, b)
    return D_ * ( b - R * parameters)
end

function solve(index)
    # A matrix of constant parameters, b is rhs vector
    A, b = gimme_gamma(index)
    parameters = zeros(Float64, MATRIX_DIMENSION)
    flag = true
    count = 0
    D_ = Diagonal(A)
    R = A - D_
    D_ = inv(D_)

    while flag
        count += 1
        parameters = jacobi(R, D_, parameters, b)
        # Quality of iteration using Frobenius norm (p=2)
        top = norm(A*parameters - b)
        bottom = norm(b)
        if top / bottom < parse(Float64, "10e-6") flag = false end
        # solution did not converge
        if isnan(top / bottom) break end
    end
    println("COUNT:", count)
    println("SOLUTION: ", A*parameters)
    println("x: ",parameters)
    println("b: ", b)
    e = abs.(A*parameters - b)
    println(round.(e; digits=7))
end

# for i in 1:3 solve(i) end

function gauss_seidel(U, L, parameters, b)
    return inv(L) * ( b - U * parameters)
end

function solve_gaus(index)
    A, b = gimme_gamma(index)
    U, L = lowerTmatrix(A)

    parameters = zeros(Float64, MATRIX_DIMENSION)
    flag = true
    count = 0
    while flag
        count += 1
        parameters = gauss_seidel(U, L, parameters, b)
        # Quality of iteration using Frobenius norm (p=2)
        top = norm(A*parameters - b)
        bottom = norm(b)
        if top / bottom < parse(Float64, "10e-6") flag = false end
        # solution did not converge
        if isnan(top / bottom) break end
    end
    println("COUNT:", count)
    println("SOLUTION: ", A*parameters)
    println("x: ",parameters)
    println("b: ", b)
    e = abs.(A*parameters - b)
    println(round.(e; digits=7))
end

solve_gaus(1)