using StatsBase

"""
    gen_moments_sel(data, mom_set, data_longmem)

Calculate moments of a time series and remove those not included in the moment 
set. If `data_longmem` is populated, long-memory moments are calculated using 
this time series instead.

# Arguments
- `data::Array{Float64,1}`: array of observations
- `mom_set::Array`: moment set
- `data_longmem=nothing`: array of observations for long-memory moments
"""
function gen_moments_sel(data::Array{Float64,1}, mom_set::Array, data_longmem=nothing)
    moments = zeros(22) # array of calculated moments

    # calculate autocorrelations of raw returns
    acf = autocor(data, 0:3)

    # calculate autocorrelations of absolute returns
    absdata = abs.(data)
    absacf = autocor(absdata, 0:101)

    # calculate autocorrelations of squared returns
    sqrdata = data.^2
    sqracf = autocor(sqrdata, 0:26)

    # FW ~ Franke & Westerhoff (2012), CL ~ Chen & Lux (2018)
    moments[1] = var(data) # variance of raw returns [CL]
    moments[2] = kurtosis(data)+3 # kurtosis of raw returns [CL]
    moments[3] = acf[2,1] # 1st lag AC of raw returns [FW,CL]
    moments[4] = acf[3,1] # 2nd lag AC of raw returns
    moments[5] = acf[4,1] # 3rd lag AC of raw returns

    moments[6] = mean(absdata) # mean of absolute returns [FW]
    moments[7] = hill(absdata, 2.5) # Hill estimator (2.5% of the right tail) of absolute returns
    moments[8] = hill(absdata, 5) # Hill estimator (5% of the right tail) of absolute returns [FW]

    moments[9] = mean(absacf[2:3,1]) # 1st lag AC of absolute returns [FW,CL]
    moments[10] = mean(absacf[5:7,1]) # 5th lag AC of absolute returns [FW,CL]
    moments[11] = mean(absacf[10:12,1]) # 10th lag AC of absolute returns [FW,CL]
    moments[12] = mean(absacf[15:17,1]) # 15th lag AC of absolute returns [CL]
    moments[13] = mean(absacf[20:22,1]) # 20th lag AC of absolute returns [CL]
    moments[14] = mean(absacf[25:27,1]) # 25th lag AC of absolute returns [FW,CL]
    moments[15] = mean(absacf[50:52,1]) # 50th lag AC of absolute returns [FW]
    moments[16] = mean(absacf[100:102,1]) # 100th lag AC of absolute returns [FW]

    moments[17] = mean(sqracf[2:3,1]) # 1st lag AC of squared returns [CL]
    moments[18] = mean(sqracf[5:7,1]) # 5th lag AC of squared returns [CL]
    moments[19] = mean(sqracf[10:12,1]) # 10th lag AC of squared returns [CL]
    moments[20] = mean(sqracf[15:17,1]) # 15th lag AC of squared returns [CL]
    moments[21] = mean(sqracf[20:22,1]) # 20th lag AC of squared returns [CL]
    moments[22] = mean(sqracf[25:27,1]) # 25th lag AC of squared returns [CL]

    # use alternative series to caculate long-memory moments
    if !isnothing(data_longmem)
        # calculate autocorrelations of absolute returns of long-memory 
        absdata_longmem = abs.(data_longmem)
        absacf_longmem = autocor(absdata_longmem, 0:101)

        # calculate autocorrelations of squared returns
        sqrdata_longmem = data_longmem.^2
        sqracf_longmem = autocor(sqrdata_longmem, 0:26)

        # FW ~ Franke & Westerhoff (2012), CL ~ Chen & Lux (2018)
        moments[11] = mean(absacf_longmem[10:12,1]) # 10th lag AC of absolute returns [FW,CL]
        moments[12] = mean(absacf_longmem[15:17,1]) # 15th lag AC of absolute returns [CL]
        moments[13] = mean(absacf_longmem[20:22,1]) # 20th lag AC of absolute returns [CL]
        moments[14] = mean(absacf_longmem[25:27,1]) # 25th lag AC of absolute returns [FW,CL]
        moments[15] = mean(absacf_longmem[50:52,1]) # 50th lag AC of absolute returns [FW]
        moments[16] = mean(absacf_longmem[100:102,1]) # 100th lag AC of absolute returns [FW]

        moments[19] = mean(sqracf_longmem[10:12,1]) # 10th lag AC of squared returns [CL]
        moments[20] = mean(sqracf_longmem[15:17,1]) # 15th lag AC of squared returns [CL]
        moments[21] = mean(sqracf_longmem[20:22,1]) # 20th lag AC of squared returns [CL]
        moments[22] = mean(sqracf_longmem[25:27,1]) # 25th lag AC of squared returns [CL]
    end

    # remove moments not included in the moment set
    moments_sel = [moments[i] for i in eachindex(mom_set) if mom_set[i] == 1]

    return moments_sel
end


"""
    gen_moments_sel_ind(data, mom_set, sample)

Calculate moments of a time series based on a sample of individual observations.

# Arguments
- `data::Array{Float64,1}`: array of observations
- `mom_set::Array`: moment set
- `sample`: sampled observations
"""
function gen_moments_sel_ind(data::Array{Float64,1}, mom_set::Array, sample)
    moments = zeros(22) # array of calculated moments

    rawdata = copy(data) # raw returns
    absdata = abs.(data) # absolute returns
    sqrdata = data.^2 # squared returns

    # sample raw returns
    rawsample = data[sample]
    rawmean = mean(rawsample)
    rawvar = var(rawsample)
    rawsampledm = rawsample .- rawmean

    # sample absolute returns
    abssample = abs.(rawsample)
    absmean = mean(abssample)
    absvar = var(abssample)
    abssampledm = abssample .- absmean

    # sample squared returns
    sqrsample = rawsample.^2
    sqrmean = mean(sqrsample)
    sqrvar = var(sqrsample)
    sqrsampledm = sqrsample .- sqrmean

    # add mean of each series to act as neutral point
    append!(rawdata, rawmean)
    append!(absdata, absmean)
    append!(sqrdata, sqrmean)

    # point samples with negative predecessors to neutral point
    sample_lag = reduce(vcat, [sample .- i for i in 1:101]')
    sample_lag[sample_lag .< 1] .= length(data) + 1

    # calculate autocorrelations of raw returns
    raw1 = rawdata[sample_lag[1, :]] .- rawmean
    raw2 = rawdata[sample_lag[2, :]] .- rawmean
    raw3 = rawdata[sample_lag[3, :]] .- rawmean

    # calculate autocorrelations of absolute returns
    abs1 = absdata[sample_lag[1:2, :]] .- absmean
    abs5 = absdata[sample_lag[4:6, :]] .- absmean
    abs10 = absdata[sample_lag[9:11, :]] .- absmean
    abs15 = absdata[sample_lag[14:16, :]] .- absmean
    abs20 = absdata[sample_lag[19:21, :]] .- absmean
    abs25 = absdata[sample_lag[24:26, :]] .- absmean
    abs50 = absdata[sample_lag[49:51, :]] .- absmean
    abs100 = absdata[sample_lag[99:101, :]] .- absmean

    # calculate autocorrelations of squared returns
    sqr1 = sqrdata[sample_lag[1:2, :]] .- sqrmean
    sqr5 = sqrdata[sample_lag[4:6, :]] .- sqrmean
    sqr10 = sqrdata[sample_lag[9:11, :]] .- sqrmean
    sqr15 = sqrdata[sample_lag[14:16, :]] .- sqrmean
    sqr20 = sqrdata[sample_lag[19:21, :]] .- sqrmean
    sqr25 = sqrdata[sample_lag[24:26, :]] .- sqrmean

    T = length(sample) # length of the sample

    # FW ~ Franke & Westerhoff (2012), CL ~ Chen & Lux (2018)
    moments[1] = var(rawsample) # variance of raw returns [CL]
    moments[2] = kurtosis(rawsample)+3 # kurtosis of raw returns [CL]
    moments[3] = ((raw1' * rawsampledm) ./ rawvar) / T # 1st lag AC of raw returns [FW,CL]
    moments[4] = ((raw2' * rawsampledm) ./ rawvar) / T # 2nd lag AC of raw returns
    moments[5] = ((raw3' * rawsampledm) ./ rawvar) / T # 3rd lag AC of raw returns

    moments[6] = mean(abssample) # mean of absolute returns [FW]
    moments[7] = hill(abssample, 2.5) # Hill estimator (2.5% of the right tail) of absolute returns
    moments[8] = hill(abssample, 5) # Hill estimator (5% of the right tail) of absolute returns [FW]
    
    moments[9] = mean((abs1 * abssampledm) ./ absvar) / T # 1st lag AC of absolute returns [FW,CL]
    moments[10] = mean((abs5 * abssampledm) ./ absvar) / T # 5th lag AC of absolute returns [FW,CL]
    moments[11] = mean((abs10 * abssampledm) ./ absvar) / T # 10th lag AC of absolute returns [FW,CL]
    moments[12] = mean((abs15 * abssampledm) ./ absvar) / T # 15th lag AC of absolute returns [CL]
    moments[13] = mean((abs20 * abssampledm) ./ absvar) / T # 20th lag AC of absolute returns [CL]
    moments[14] = mean((abs25 * abssampledm) ./ absvar) / T # 25th lag AC of absolute returns [FW,CL]
    moments[15] = mean((abs50 * abssampledm) ./ absvar) / T # 50th lag AC of absolute returns [FW]
    moments[16] = mean((abs100 * abssampledm) ./ absvar) / T # 100th lag AC of absolute returns [FW]

    moments[17] = mean((sqr1 * sqrsampledm) ./ sqrvar) / T # 1st lag AC of squared returns [CL]
    moments[18] = mean((sqr5 * sqrsampledm) ./ sqrvar) / T # 5th lag AC of squared returns [CL]
    moments[19] = mean((sqr10 * sqrsampledm) ./ sqrvar) / T # 10th lag AC of squared returns [CL]
    moments[20] = mean((sqr15 * sqrsampledm) ./ sqrvar) / T # 15th lag AC of squared returns [CL]
    moments[21] = mean((sqr20 * sqrsampledm) ./ sqrvar) / T # 20th lag AC of squared returns [CL]
    moments[22] = mean((sqr25 * sqrsampledm) ./ sqrvar) / T # 25th lag AC of squared returns [CL]

    # remove moments not included in the moment set
    moments_sel = [moments[i] for i in eachindex(mom_set) if mom_set[i] == 1]

    return moments_sel
end

"""
    hill(data, pct)

Calculate Hill estimator at `pct` of the right tail.

# Arguments
- `data::Array{Float64,1}`: array of observations
- `pct`: percentage of the right tail
"""
function hill(data::Array{Float64,1}, pct)
    k = floor(Int, length(data)/100*pct) # determine number of considered points

    sorted = sort(data, rev=true) # sort data from highest to lowest values
    res = sorted[1:k]/sorted[k+1] # normalize considered points

    return ((1/k)*sum(log.(res)))^-1
end
