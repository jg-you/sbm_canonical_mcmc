# sbm_canonical_mcmc

C++ implementation of a MCMC sampler for the (canonical) MCMC.

## Usage

Depends on `boost::program_options` and `cmake`.

Compilation:

	cmake .
	make

The binaries are built in `bin/`.

Example call:

	bin/mcmc -e example_edge_list.txt -P 0.6 0.1 0.1 0.6 -n 20 20 -r -t 100 -f 40 -b 50 -s

For the edge list at `example_edge_list.txt`, sample from the marginal of the SBM with probability matrix (`-P`)

	0.6  0.1
	0.1  0.6

and two blocks of 20 nodes (`-n 20 20`), with randomized initial condition (`-r`), every 100 MCMC moves (`-t 100`), with a 
burn-in period  of 50 steps (`-b 50`), every 40 moves (`-f 40`).
Uses the single vertex move proposal distribution (`-s`).

If the build was succesfull, the output should look like

    edge_list_path: example_edge_list.txt
    probabilities:
    0.6 0.1 
    0.1 0.6 
    sizes (g=2): 20 20 
    burn_in: 50
    sampling_steps: 100
    sampling_frequency: 40
    maximize: false
    use_ppm: false
    use_single_vertex: true
    randomize: true
    seed: 1434121506
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    acceptance ratio 0.5055


## Companion article

Please cite

**Finite size analysis of the detectability limit of the stochastic block model**

Jean-Gabriel Young, Patrick Desrosiers, Laurent Hébert-Dufresne, Edward Laurence, Louis J. Dubé

[arXiv:1701.00062](https://arxiv.org/abs/1701.00062)
