# Moment set selection for the SMM using simple machine learning

This repository contains the official implementation of [Moment set selection for the SMM using simple machine learning](https://doi.org/10.1016/j.jebo.2023.05.040) by [Eric Zila](zila.eric@gmail.com) and [Jiri Kukacka](jiri.kukacka@fsv.cuni.cz).

If you find our work useful, we encourage you to use the following citation:
```
@article{zilakukacka2023,
    title = {Moment set selection for the {SMM} using simple machine learning},
    author = {Eric Zila and Jiri Kukacka},
    journal = {Journal of Economic Behavior & Organization},
    volume = {212},
    pages = {366-391},
    year = {2023},
    issn = {0167-2681},
    doi = {https://doi.org/10.1016/j.jebo.2023.05.040},
    url = {https://www.sciencedirect.com/science/article/pii/S0167268123001944}
}
```

## Requirements

To download the content of the repository and move to the root directory, run:
```
git clone git@github.com:jirikukacka/Zila_Kukacka_Moment_Selection.git
cd Zila_Kukacka_Moment_Selection/
```

To install the necessary packages, run [Julia]((https://julialang.org/)) (>= 1.6.1) in the root directory and execute:
```
julia> using Pkg
julia> Pkg.activate(pwd())
julia> Pkg.instantiate()
```

## Reproduction

The repository contains all the code needed to reproduce experiments from the paper for the Markov-switching multifractal (MSM) model. The necessary scripts are located in the `scripts/` folder.

To estimate the MSM model using the three benchmark moment sets employed in the paper, execute:
```
julia scripts/main_markov256state_bench.jl
```

To perform the backward stepwise moment elimination (BSME) algorithm as proposed in the paper to the MSM model, execute:
```
julia scripts/main_markov256state_bsw.jl
```

To perform the forward stepwise moment selection (FSMS) algorithm as proposed in the paper to the MSM model, execute:
```
julia scripts/main_markov256state_bsw.jl
```

## Results

For the Markov-switching multifractal model, the full set of results is located in the `results/` folder. The corresponding visualizations are located in the `plots/` folder. In order to reproduce them, please, refer to the `notebooks/markov256state.ipynb` Jupyter notebook file executable with a Julia kernel.

## Data

Empirical data used for empirical estimation are located in the `data/` folder. 

## License

All content in the repository is licensed under the MIT license. More information can be found in the [license file](LICENSE).
