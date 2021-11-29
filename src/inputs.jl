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
    if !isa(index, Integer) || index < 1 || index > 3
        throw(DomainError(index, "Gamma can be only 5, 2, 0.5 with respective indices 1, 2, 3."))
    end
    gamma = [5.0, 2.0, 0.5]

    dv = fill(gamma[index], MATRIX_DIMENSION)
    ev = fill(-1.0, MATRIX_DIMENSION - 1)
    Au = Bidiagonal(dv, ev, :U)
    Av = Bidiagonal(zeros(Float64, MATRIX_DIMENSION), ev, :L)
    A = Au + Av

    b = zeros(Float64, MATRIX_DIMENSION)
    for i in 1:MATRIX_DIMENSION
        b[i] = gamma[index] - 2
    end
    # first and last element in b vector is only gamma - 1
    b[1] = gamma[index] - 1
    b[MATRIX_DIMENSION] = gamma[index] - 1
    return A, b
end

