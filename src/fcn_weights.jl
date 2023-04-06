using LinearAlgebra
using Distributions
using Random
using StatsBase

include("fcn_moments.jl")

"""
    gen_weights(setup, data, mom_set)

Generate SMM weighting matrix.

# Arguments
- `setup::Dict`: full setup dictionary
- `data::Array{Float64}`: (pseudo-)empirical dataset
- `mom_set::Array`: moment set
"""
function gen_weights(setup::Dict, data::Array{Float64}, mom_set::Array)
    # identity matrix
    if setup["wgt"]["method"] == "eye"
        return weights_identity(mom_set)
    # non-overlapping block bootstrap matrix
    elseif setup["wgt"]["method"] == "fw2012"
        return weights_fw2012(data, setup["wgt"]["bootsize"], mom_set)
    # overlapping block bootstrap matrix
    elseif setup["wgt"]["method"] == "fw2012overlap"
        return weights_fw2012(data, setup["wgt"]["bootsize"], mom_set, overlap=true)
    # diagonal overlapping block bootstrap matrix
    elseif setup["wgt"]["method"] == "fw2012overlapdiag"
        return weights_fw2012(data, setup["wgt"]["bootsize"], mom_set, overlap=true, diag=true)
    # history sampling bootstrap matrix
    elseif setup["wgt"]["method"] == "fw2016"
        return weights_fw2016(data, setup["wgt"]["bootsize"], mom_set)
    else
        error("The chosen method of computing weighting matrix is unknown.")
    end
end


"""
    weights_identity(mom_set)

Generate identity weighting matrix.

# Arguments
- `mom_set::Array`: moment set
"""
weights_identity(mom_set::Array) = Matrix{Float64}(I, sum(mom_set), sum(mom_set))


"""
    weights_fw2012(data, bootsize, mom_set[, overlap=false, diag=false])

Generate (non-)overlapping block bootstrap weighting matrix.

Implementation details can be found in Franke & Westerhoff (2012). Blocks of 250
observations are used for short-memory moments, while blocks of 750 observations
are used for long-memory moments.
    
If `overlap=false`, the original version with non-overlapping blocks is computed
in line with Franke & Westerhoff (2012). If `overlap=true`, an adjusted version
with overlapping blocks is computed according to Zila & Kukacka (2023).

If `diag=false`, the full weighting matrix is returned. If `diag=true`,
non-diagonal elements of the matrix are removed.

# Arguments
- `data::Array{Float64}`: (pseudo-)empirical dataset
- `bootsize::Int`: bootstrap repetitions
- `mom_set::Array`: moment set
- `overlap::Bool=false`: flag on whether blocks should overlap or not
- `diag::Bool=false`: flag on whether matrix should be diagonalized or not
"""
function weights_fw2012(data::Array{Float64}, bootsize::Int, mom_set::Array; overlap::Bool=false, diag::Bool=false)
    moments_bt = zeros(sum(mom_set), bootsize) # array of bootstrapped moments

    # calculate number of non-overlapping blocks per boostrapped samples
    blockcount_250 = floor(Int, length(data)/250) # blocks of 250 observations
    blockcount_750 = floor(Int, length(data)/750) # blocks of 750 observations

    for i in 1:bootsize
        # draw starting indices of overlapping blocks
        if overlap
            init_250 = rand(DiscreteUniform(1, size(data)[1]-249), blockcount_250)
            init_750 = rand(DiscreteUniform(1, size(data)[1]-749), blockcount_750)
        # draw starting indices of non-overlapping blocks
        else
            init_250 = rand(DiscreteUniform(0, blockcount_250-1), blockcount_250)
            init_750 = rand(DiscreteUniform(0, blockcount_750-1), blockcount_750)
        end

        # prepare arrays for bootstrapped data
        cur_data_250 = Float64[]
        cur_data_750 = Float64[]

        # fill array with blocks of 250 observations
        for j in init_250
            if overlap
                append!(cur_data_250, data[j:j+249])
            else
                append!(cur_data_250, data[j*250+1:j*250+250])
            end
        end

        # fill array with blocks of 750 observations
        for j in init_750
            if overlap
                append!(cur_data_750, data[j:j+749])
            else
                append!(cur_data_750, data[j*750+1:j*750+750])
            end
        end

        # calculate moments from the bootstrapped series
        moments_bt[:,i] = gen_moments_sel(cur_data_250, mom_set, cur_data_750)
    end

    weights = inv(cov(transpose(moments_bt))) # calculate weights

    if diag
        weights = diagm(diag(weights)) # remove non-diagonal elements
    end

    return weights
end


"""
    weights_fw2016(data, bootsize, mom_set[, diag=false])

Generate history sampling bootstrap weighting matrix.

Implementation details can be found in Franke & Westerhoff (2016).

If `diag=false`, the full weighting matrix is returned. If `diag=true`,
non-diagonal elements of the matrix are removed.

# Arguments
- `data::Array{Float64}`: (pseudo-)empirical dataset
- `bootsize::Int`: bootstrap repetitions
- `mom_set::Array`: moment set
- `diag::Bool=false`: flag on whether matrix should be diagonalized or not
"""
function weights_fw2016(data::Array{Float64}, bootsize::Int, mom_set::Array, diag::Bool=false)
    moments_bt = zeros(sum(mom_set), bootsize) # array of bootstrapped moments

    # draw indices of sampled observations for all repetitions
    samples = rand(DiscreteUniform(1, length(data)), bootsize, length(data))

    # calculate moments from each repetition's sample
    for i in 1:bootsize
        moments_bt[:, i] = gen_moments_sel_ind(data, mom_set, samples[i, :])
    end

    moments_bt_demean = moments_bt .- mean(moments_bt, dims=2) # demean moments

    weights = inv(cov(moments_bt_demean, dims=2)) # calculate weights

    if diag
        weights = diagm(diag(weights)) # remove non-diagonal elements
    end

    return weights
end
