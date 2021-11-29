#=
main:
- Julia version: 
- Author: sutymate
- Date: 2021-11-27
=#

include("../src/inputs.jl")


# first try, not working
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


function gauss_seidel(U, L_, parameters, b)
    return L_ * ( b - U * parameters)
end

function make_graphs(jc_e, gs_e, index)
    gamma = [5.0, 2.0, 0.5]
    # returned tuple from solve
    jc_e = jc_e[3]
    gs_e = gs_e[4]
    println(jc_e)
    println(gs_e)
    # number of iterations
    its = length(gs_e)
    if length(jc_e) > length(gs_e) its = length(jc_e) end
    p = plot([
        scatter(x=1:its, y=jc_e, mode="lines", name="Jacobi", color=:red),
        scatter(x=1:its, y=gs_e, mode="marker", name="Gauss-Seidel", color=:blue)],
        Layout(title="Norm error for gamma "*string(gamma[index]), yaxis_title="Norm error", xaxis_title="Number of iterations"))
    savefig(p, string(index)*"_.png")

end

function solve(index, calculation_method="gs", verbose=true)
    jc = false
    gs = false
    if calculation_method == "jacobi" || calculation_method == "gs"
        if calculation_method == "jacobi"
             jc = true
        else gs = true end
    else throw(ArgumentError("Calculation method can be one of: \"jacobi\", \"gs\".")) end

    # A matrix of constant parameters, b is vector from A*x=b
    A, b = gimme_gamma(index)
    parameters = zeros(Float64, MATRIX_DIMENSION)
    flag = true
    iterations = 0

    if jc
        method_name = "Jacobi"
        D_ = Diagonal(A)
        R = A - D_
        D_ = inv(D_)
    end
    if gs
        method_name = "Gauss-Seidel"
        A, b = gimme_gamma(index)
        U, L = lowerTmatrix(A)
        L_ = inv(L)
    end

    if verbose println("Calculating using "*method_name*".") end

    top = 0
    bottom = 0
    gs_e = []
    jc_e = []
    # set upper limit of iterations to 10 000
    while flag && iterations < 10^4
        iterations += 1
        if jc parameters = jacobi(R, D_, parameters, b) end
        if gs parameters = gauss_seidel(U, L_, parameters, b) end

        # Quality of iteration using Frobenius norm (p=2)
        top = norm(A*parameters - b)
        bottom = norm(b)
        if top / bottom < parse(Float64, "10e-6") flag = false end
        # solution did not converge
        if verbose
            e = abs.(A*parameters - b)
            println("Iteration: "*string(iterations)*", error: ", round(top / bottom; digits=7))
        end
        if jc push!(jc_e, top/bottom)
        else push!(gs_e, top/bottom) end
        if isnan(top / bottom) break end
    end
    norm_e = top / bottom
    if isnan(norm_e)
        println("\n\tCALCULATION DID NOT CONVERGE after "*string(iterations)*" iterations.\n")
        return -1, -1, -1, -1
    end
    if verbose
        println("Iterations:", iterations)
        println("SOLUTION x: ",parameters)
        println("CALCULATED RESULT A*x: ", round.(A*parameters; digits=3))
        println("GIVEN RESULT b: ", b)
        e = abs.(A*parameters - b)
        println("Calculated error: ", round.(e; digits=7))
    end
    if ! verbose println("CALCULATED RESULT A*x: ", round.(A*parameters; digits=3)) end
    println("Normative error: ", norm_e)
    println()
    return norm_e, iterations, jc_e, gs_e
end

# select one method
calculation_method = ["jacobi", "gs"]
cm = "jacobi"
# select gamma index
gamma_index = [1, 2, 3]
index = 1
# choose if you want more detailed output
verbose = true
# norm_e is final value of criterial function (eucl. norm
# ||Ax_ - b|| / || b ||
# iterations is number of iterations
# jc_e, resp. gs_e, is a list of criterial function values for each iteration
norm_e, iterations, jc_e, gs_e = solve(index, cm, verbose)



### REMOVE comment for making graphs
# using PlotlyJS
# e = Dict()
# for cm in calculation_method
#     e[cm] = Dict()
#     for i in 1:3
#         e[cm][i] = solve(i, cm, false)
#     end
# end

# for i in 1:3
#     make_graphs(e["jacobi"][i], e["gs"][i], i)
# end

# for method in keys(e)
#     println(method)
#     for idx in keys(e[method])
#         println(idx, "\t", e[method][idx])
#     end
# end


