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

function solve()
    # A matrix of constant parameters, b is rhs vector
    A, b = gimme_gamma(1)
    prev_parameters = zeros(Float64, MATRIX_DIMENSION)
    next_parameters = zeros(Float64, MATRIX_DIMENSION)
    flag = true
    count = 0
    println(A)
    println(b)
    while flag
        count += 1
        prev_parameters = next_parameters
        for index in 1:MATRIX_DIMENSION
            row = A[index, :]
            # calculate x_index
            next_parameters[index] = apx(index, row, prev_parameters, b)
        end
        # Quality of iteration using Frobenius norm (p=2)
        top = norm(A*next_parameters - b)
        bottom = norm(b)
        if top / bottom < parse(Float64, "10e-6") flag = false end
        # solution did not converge
        if isnan(top / bottom) break end
    end
    println("COUNT:", count)
    println("SOLUTION: ", A*next_parameters)
    println("x: ",next_parameters)
    println("b: ", b)
end

solve()