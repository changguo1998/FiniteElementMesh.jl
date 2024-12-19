function _polycoef_1d_linear(order::Int)
    if order <= 1
        error("order <= 1")
    end
    return collect(range(-1.0, 1.0, length=order+1))
end

function _polycoef_1d_gll end

function _polycoef_2d_combine_triangle(f_coef::Function, order::Int)
    coef1d = f_coef(order)
    np = count(c -> (coef1d[c.I[1]] + coef1d[c.I[2]]) <= 1.0, CartesianIndices((length(coef1d), length(coef1d))))
    coef2d = zeros(2, np)
    p = 1
    for i = eachindex(coef1d), j = eachindex(coef1d)
        if (coef1d[i] + coef1d[j]) > 1.0
            continue
        end
        coef2d[1, p] = coef1d[i]
        coef2d[2, p] = coef1d[j]
        p += 1
    end
    return coef2d
end

function _polycoef_2d_combine_rectangle(f_coef::Function, order::Int)
    coef1d = f_coef(order)
    coef2d = zeros(2, length(coef1d) * length(coef1d))
    p = 1
    for i = eachindex(coef1d), j = eachindex(coef1d)
        coef2d[1, p] = coef1d[i]
        coef2d[2, p] = coef1d[j]
        p += 1
    end
    return coef2d
end
