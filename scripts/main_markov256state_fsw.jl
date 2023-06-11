# paralell computing
using Distributed
addprocs(48)

# initialize project using DrWatson
@everywhere using DrWatson
@everywhere @quickactivate "Zila_Kukacka_2023"

include(srcdir("fcn_smm_ml.jl"))

# machine learning setup ~ {options}
ml_set = Dict("method" => "fsw", # machine learning method ~ {"fsw","bsw","sms"}
              "bench" => nothing, # benchmark set if "method" == "sms" ~ {nothing,"chl4","chl15","fw9","full"}
              "set" => nothing) # moment set if "bench" == nothing

# simulated method of moments setup [default]
smm_set = Dict("rep" => 96, # number of repetitions [96]
               "emp" => nothing, # source of empirical data ["SP500_1980_2022_logret.jld"]
               "simfactor" => 1) # simulated series length factor [1]

# model setup [default] ~ {options}
mod_set = Dict("model" => "markov256state", # model name
               "obs" => 6750, # number of observations [6750]
               "burn" => 200, # burn-in period length [200]
               "cali" => "lux2022", # model calibration
               "cons" => "lux2022") # search constraints

# optimisation setup [default]
opt_set = Dict("inits" => 1, # number of initial points [1]
               "sim" => 100, # number of simulations [100]
               "iter" => 2000) # number of iterations [4000]

# weighting matrix setup [default] ~ {options}
wgt_set = Dict("method" => "fw2012overlap", # approach to weighting matrix construction ~ {"eye","fw2012","fw2012overlap","fw2012overlapdiag","fw2016"}
               "blocksize" => 250, # block size for block-bootstrapped weighting matrix [250]
               "bootsize" => 5000, # repetitions for block-bootstrapped weighting matrix [5000]
               "blockcount" => 1000) # number of pasted blocks for long block-bootstrapped weighting matrix [1000]

# full setup dictionary
setup = Dict("ml" => ml_set, # machine learning setup
             "smm" => smm_set, # simulated method of moments setup
             "mod" => mod_set, # model setup
             "opt" => opt_set, # optimisation setup
             "wgt" => wgt_set) # weighting matrix setup

@time begin
    results = smm_init(setup) # estimation initialization
end
