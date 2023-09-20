function update_bias_parametric!(b::Metadynamics, cv, itrj)
    batchsize = b.batchsize
    push!(b.cv_storage, cv)
    push!(b.bias_storage, b(cv))
    push!(b.fullbias_storage, b.values)
    update_bias_regular!(b, cv)

    if b.symmetric
        update_bias_regular!(b, -cv)
    end

    if length(b.cv_storage)%batchsize==0
        p = b.current_parameters
        f = OptimizationFunction(b.testfun, Optimization.AutoForwardDiff())
        sol = solve(
            OptimizationProblem(f, p, b, lb=b.lower_bounds, ub=b.upper_bounds),
            b.minimizer,
            time_limit = 60,
        )
        b.current_parameters .= sol.u
        test = b.testfun(p, b)
        println_verbose(
            b.KS_fp,
            itrj, " ", sol.u[1], " ", sol.u[2], " ", sol.u[3], " ", test,
            "# itrj parameters $(b.testfun)"
        )
        println("==============================================")
        flush(b.KS_fp)
    end

    return nothing
end

function parametricFES(parameters, cv)
    A, B, C = parameters
    return -A * cv^2 - B * cos(π * cv * C)^2
end

function parametricBias(parameters, cv)
    A, B, C = parameters
    return A * cv^2 + B * cos(π * cv * C)^2
end

function KStest(parameters, b::Metadynamics)
    #Sg = cv -> parametricFES(parameters, cv)
    batchsize = b.batchsize

    cv_data = b.cv_storage
    len = length(cv_data)
    itvl = Int(len / batchsize)

    weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
    bias_data = b.bias_storage[1:itvl:end]
    fullbias_data = b.fullbias_storage[1:itvl:len]
    cv_data = cv_data[1:itvl:end]
    cdf = ecdf(cv_data, weights=weight.(cv_data, bias_data, fullbias_data))
    #cdf = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
    sorteddata = cdf.sorted_values

    cvmin = sorteddata[1]
    cvmax = sorteddata[end]

    l = abs(0.0 - (sorteddata[1] - cvmin) / (cvmax - cvmin))
    r = abs(cdf(sorteddata[1]) - (sorteddata[1] - cvmin) / (cvmax - cvmin))

    for i in 2:batchsize
        P_i = (sorteddata[i] - cvmin) / (cvmax - cvmin)
        lnext = abs(cdf(sorteddata[i-1]) - P_i)
        rnext = abs(cdf(sorteddata[i]) - P_i)
        l = l ≤ lnext ? lnext : l
        r = r ≤ rnext ? rnext : r
    end

    return max(l, r)
end

function GADtest(parameters, b::Metadynamics)
    #Sg = cv -> parametricFES(parameters, cv)
    batchsize = b.batchsize

    cv_data = b.cv_storage
    len = length(cv_data)
    itvl = Int(len / batchsize)

    weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
    bias_data = b.bias_storage[1:itvl:len]
    fullbias_data = b.fullbias_storage[1:itvl:len]
    cv_data = cv_data[1:itvl:end]

    #F = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
    F = ecdf(cv_data, weights = weight.(cv_data, bias_data, fullbias_data))
    sorteddata = F.sorted_values
    #cvmin, cvmax = get_CVlims(b)
    cvmin = floor(minimum(sorteddata))
    cvmax = ceil(maximum(sorteddata))

    S = 0.0

    for i in 1:batchsize - 1
        p_i = F(sorteddata[i])
        u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
        u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
        S += p_i^2 * log(u_ip1 / u_i) - (p_i-1)^2 * log((1-u_ip1) / (1-u_i))
    end

    u_1 = wantedCDF(sorteddata[1], cvmin, cvmax)
    u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
    return batchsize * (-1 - log(u_n) - log(1-u_1) + S)
end

function GADLTtest(parameters, b::Metadynamics)
    #Sg = cv -> parametricFES(parameters, cv)
    batchsize = b.batchsize

    cv_data = b.cv_storage
    len = length(cv_data)
    itvl = Int(len / batchsize)

    weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
    bias_data = b.bias_storage[1:itvl:len]
    fullbias_data = b.fullbias_storage[1:itvl:len]
    cv_data = cv_data[1:itvl:end]

    #F = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
    F = ecdf(cv_data, weights = weight.(cv_data, bias_data, fullbias_data))
    sorteddata = F.sorted_values
    cvmin = floor(minimum(sorteddata))
    cvmax = ceil(maximum(sorteddata))

    S = 0.0

    for i in 1:batchsize - 1
        p_i = F(sorteddata[i])
        u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
        u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
        S += p_i^2 * log(u_ip1/u_i) + 2p_i * (u_i - u_ip1)
    end

    p_n = F(sorteddata[end])
    u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
    return batchsize * (0.5 - 2p_n * (1-u_n) - p_n^2 * log(u_n) + S)
end

function weight_for_cdf(parameters, cv, bias, fullbias, b::Metadynamics)
    cv_vals = b.cv_vals
    inorm = parametric_norm(parameters, fullbias, cv_vals)
    weight = inorm * exp(parametricFES(parameters, cv) + bias)
    return weight
end

function parametric_norm(parameters, old_bias, cv_vals)
    normA = 0.0
    i = 0

    for cv in cv_vals
        i += 1
        normA += exp(-parametricFES(parameters, cv) - old_bias[i])
        #normA += exp( -parametricFES(parameters, cv) -
        #	parametricBias(old_bias, cv) )
    end

    return normA
end

function wantedCDF(cv, cvmin, cvmax)
    return (cv - cvmin) / (cvmax - cvmin)
end

function parametric_to_bias!(b::Metadynamics)
    parameters = b.current_parameters

    for (idx, cv) in enumerate(b.cv_vals)
        b[idx] = parametricBias(parameters, cv)
    end

    return nothing
end
