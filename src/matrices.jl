"""
$(TYPEDSIGNATURES)

Outputs the classic standard mass isotopomer distribution (MID) for a fragment 
    with N atoms of the given element, using combinatorial probabilities.
Not to be used for atoms labeled through an isotope tracer experiment, e.g. 
    backbone carbon atoms in a 13C tracer experiment.

# Arguments
- `isotopes::Dict`: dictionary containing the natural abundances of the relevant isotopes.
- `element`: the element the standard MID should be computed for.
- `N::Int`: the number atoms of the given element.
"""
function standard_element_MID(isotopes::Dict, element, N::Int)
    # function to calculate a multinomial coefficient, n should be a number and k a vector
    multinomial_coef(n, k) = factorial(Int(n)) / cumprod(factorial.(Int.(k)))[end]

    result = [] # initialize results object

    # for H, N or unlabeled C
    if element == "H" || element == "N" || element == "C"
        result = [binomial(N, i) * (isotopes[element][1]^(N - i) * isotopes[element][2]^i) for i = 0:(N)]

    else
        for i = 0:N
            push!(result, 0)
            for j = 0:N, k = 0:N
                # for O & Si
                if (element == "O" || element == "Si") && k >= 0 && i == j + 2 * k
                    c = multinomial_coef(N, [j, k, N - j - k])
                    ip = isotopes[element][1]^(N - j - k) *
                         isotopes[element][2]^j *
                         isotopes[element][3]^k
                    result[i+1] += c * ip
                else
                    for z = 0:N
                        # for S
                        if element == "S" && k >= 0 && i == j + 2 * k + 4 * z
                            c = multinomial_coef(N, [j, k, z, N - j - k - z])
                            ip = isotopes[element][1]^(N - j - k - z) *
                                 isotopes[element][2]^j *
                                 isotopes[element][3]^k *
                                 isotopes[element][4]^z
                            result[i+1] += c * ip
                        end
                    end
                end
            end
        end
    end
    return result
end


"""
$(TYPEDSIGNATURES)

Outputs the correction matrix for a fragment with N atoms of the given element.
All correction matricies will be l+1 x l+1 with l being the total number of possibly labeled atoms.

# Arguments
- `element`: the element the standard MID should be computed for.
- `N::Int`: the number atoms of the given element.
- `l::Int`: total number of possibly labeled atoms, default is -1, only needed for unlabeled atoms.
- `tracer_purity`: purity of the tracer, default is 1 (100%).
"""
function element_CM(element, N::Int; l::Int = -1, tracer_purity = 1.0, isotope_Dict = isotope_NAs)
    # TODO: change l to not keyword arg and test if N == l for labeled C
    # dictionary of common natural isotope abundances

    # labeled C
    if element == "LabC"
        CM = Matrix{Float64}(undef, N + 1, 0)
        for n = 0:N
            MID = [binomial(N - n, i) * (isotope_Dict["C"][1]^(N - n - i) * isotope_Dict["C"][2]^i) for i = 0:(N-n)]
            n_result = vcat(zeros(Int64, 1, n)', MID)

            # tracer purity
            tracer = [1.0 - tracer_purity, tracer_purity]
            n_result = conv(n_result, tracer)

            CM = hcat(CM, n_result[2:end])
        end

    # unlabeled C and all other elements
    elseif l >= 0
        CM = Matrix{Float64}(undef, l + 1, 0)
        standard_0 = vcat(standard_element_MID(isotope_Dict, element, N), zeros(Int64, 1, l)')
        for j = 0:l
            CM = hcat(CM, circshift(standard_0, j)[1:l+1])
        end
    end
    return Float64.(CM)
end

"""
$(TYPEDSIGNATURES)

Outputs the full correction matrix for a fragment.
The correction matricies will be l+1 x l+1 with l being the total number of possibly labeled atoms.

# Arguments
- `formula`: chemical formula of the fragment to be corrected.
- `tracer_purity`: purity of the tracer, default is 1 (100%).
"""
function fragment_CM(formula; tracer_purity = 1.0)
    # create a dictionary from the chemical formula of the element
    elements = split(formula, r"[0-9]", keepempty = false)
    amounts = split(formula, r"[aA-zZ]", keepempty = false)
    fragment_dict = Dict(elements .=> map(x -> parse(Int, x), amounts))
    fragment_dict["C"] = fragment_dict["C"] - fragment_dict["LabC"]

    # multiplies the individual element CMs to create the full CM
    full_CM = 1
    for (key, item) in fragment_dict
        full_CM *= element_CM(key, item, l = fragment_dict["LabC"], tracer_purity = tracer_purity)
    end
    return full_CM
end

