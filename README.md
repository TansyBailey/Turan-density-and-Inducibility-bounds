# Density functions

Required to run `density_bounds.sage` and `density_bounds_induced.sage`:
- SageMath (used SageMath version 10.0)
- CVXOPT

## To run the code

For the main script `density_bounds.sage`, run 

`main_function(r,k,F,P,problem)`

where
- `r` is a scalar (as in $r$-uniform hypergraph),
- `k` is a scalar (as in number of vertices),
- `F` is a list of incidence structures (list of forbidden graphs $\mathcal{F}$),
- `P` is a list of incidence structures (list of graphs in inducibility problem),
- `problem` is a string: `'D'` or `'I'`, indicating whether the problem is of Turán density in nature or of inducibility in nature, respectively.

See `example.txt` for examples on this code execution.

## Purpose of the code

The code in the file `density_bounds.sage` is written to produce bounds on inducibility and Turán density of given graphs/graph families.


The main function relies on 3 preliminary functions to compute the output: `F_free_rgraphs`, `flag_isomorphic`, and `flag_density`.

## Input information

- `F_free_rgraphs(r,k,F)`
    - `r`: a scalar (as in an $r$-uniform hypergraph),
    - `k`: a scalar (number of vertices),
    - `F`: a list of incidence structures (list of forbidden graphs $\mathcal{F}$).
- `flag_isomorphic(F,G)`
    - `F`: an incidence structure (a flag),
    - `G`: an incidence structure (a flag).
- `flag_density(F,G,H)`
    - `F`: an incidence structure (a flag),
    - `G`: an incidence structure (a flag),
    - `H`: an incidence structure (a graph).



## Output information

- `F_free_rgraphs(r,k,F)`
    - Outputs a list of incidence structures (list of $\mathcal{F}$-free graphs).
- `flag_isomorphic(F,G)`
    - Outputs a boolean determining whether $F$ and $G$ are isomorphic as flags.
- `flag_density(F,G,H)`
    - Outputs a scalar value $p(F,G;H)$.


## Storage of outputs

Overall, the code produces 4 text files containing information about the input problem:

- `F-free graphs.txt`
- `flags.txt`
- `probability matrices.txt`
- `output.txt`

where `output.txt` contains the information regarding the number of graphs, types and flags, and the SDP and its output. 

The user must provide the name for `savedirectory` in order to determine where these files are saved. If they do not wish to save these files, their commands can be commented out in `main_function`.

## Secondary function

There is also the function `density_bounds_induced.sage`. The inputs of all functions in this file are identical to that of the main code. However, we have a new subfunction `induced_F_free_rgraphs` that generates all induced $\mathcal{F}$-free graphs. The main function of this code, `induced_main_function`, is executed in the same way as `main_function`, but uses `induced_F_free_rgraphs` at points in its body, and minimises the objective of the inducibility SDP. 
