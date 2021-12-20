#=
inputs:
- Julia version: 
- Author: sutymate
- Date: 2021-11-27
=#

using LinearAlgebra
const  MATRIX_DIMENSION = 20


function gimme_gamma(index)
    #=
    index is index of ordered set [5, 2, 0.5]
    gamma values.
    Creates a 20x20 matrix with double precision Float64,
    where all values are -1, except for diagonal. Diagonal is
    one from [5, 2, 0.5] defined by index in parameter.
    =#
#     for i in 1:10 println(parse(Float64, "10e-2")) end
    if !isa(index, Integer) || index < 1 || index > 3
        throw(DomainError(index, "Gamma can be only 5, 2, 0.5 with respective indices 1, 2, 3."))
    end
    gamma = [5.0, 2.0, 0.5]

    # gamma value on diagonal
    dv = fill(gamma[index], MATRIX_DIMENSION)
    # -1 on sides of diagonal
    ev = fill(-1.0, MATRIX_DIMENSION - 1)
    # other positions are zero
    A = SymTridiagonal(dv, ev)
    # b is filled with gamma - 2
    b = fill(gamma[index] - 2, MATRIX_DIMENSION)
    # first and last element in b vector is gamma - 1
    b[1] = gamma[index] - 1
    b[MATRIX_DIMENSION] = gamma[index] - 1
    return A, b
end


function lowerTmatrix(A)
    n, m = size(A)
    if n != m throw(DomainError(index, "Matrix is not symmetric.")) end
    # prepare Upper Triangular and Lower Triangula matrices
    U = zeros(Float64, (n, n))
    L = zeros(Float64, (n, n))
    for row in 1:n
        for column in 1:n
            if row >= column
                L[row, column] = A[row, column]
            else
                U[row, column] = A[row, column]
            end
        end
    end
    # A = U + L
    return U, L
end
