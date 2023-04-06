using LaTeXStrings

include("fcn_watson.jl")
include(modelsdir("markov.jl"))

# named moment sets
SET = Dict(
    # Chen & Lux (2018) ~ 4 moments
    "chl4" =>  [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0], 
    # Chen & Lux (2018) ~ 15 moments
    "chl15" => [1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1], 
    # Franke & Westerhoff (2012) ~ 9 moments
    "fw9" =>   [0,0,1,0,0,1,0,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0], 
    # full moment set ~ 22 moments
	"full" =>  [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], 
)

# Markov-switching model
MARKOV = Dict(
    # calibrations
    "cali" => Dict(
        "lux2022" => [1.2, 1.0], # Lux (2022)
    ),
    # search constraints
    "cons" => Dict(
        "lux2022" => [(1.0, 2.0),  (0.0, 5.0)], # Lux (2022)
    ),
    "gnam" => [L"m_0", L"\sigma"], # parameter names
    "gdim" => (1, 2), # parameter plot dimensions
    "gord" => [1, 2], # parameter plot order
)

# implemented models' parameter settings
MODICT = Dict(
    "markov4state" => MARKOV,
    "markov16state" => MARKOV,
    "markov256state" => MARKOV,
)

# implemented models
IMPDICT = Dict(
    "markov4state" => markov4state,
    "markov16state" => markov16state,
    "markov256state" => markov256state,
)


"""
    get_model_func(mod_set)

Return model's implementation based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_func(mod_set::Dict) = IMPDICT[mod_set["model"]]


"""
    get_model_cali(mod_set)

Return model's calibration based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_cali(mod_set::Dict) = MODICT[mod_set["model"]]["cali"][mod_set["cali"]]


"""
    get_model_cons(mod_set)

Return model's search constraints based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_cons(mod_set::Dict) = MODICT[mod_set["model"]]["cons"][mod_set["cons"]]


"""
    get_model_gnam(mod_set)

Return model's parameter names based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_gnam(mod_set::Dict) = MODICT[mod_set["model"]]["gnam"]


"""
    get_model_gdim(mod_set)

Return model's parameter plot dimensions based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_gdim(mod_set::Dict) = MODICT[mod_set["model"]]["gdim"]


"""
    get_model_gord(mod_set)

Return model's parameter plot order based on model setup dictionary.

# Arguments
- `mod_set`: model setup dictionary
"""
get_model_gord(mod_set::Dict) = MODICT[mod_set["model"]]["gord"]
