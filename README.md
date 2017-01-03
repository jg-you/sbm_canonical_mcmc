# sbm_canonical_mcmc

C++ implementation of a MCMC sampler for the (canonical) MCMC.

## Usage

### Compilation

Depends on `boost::program_options` and `cmake`.

Compilation:

	cmake .
	make

The binaries are built in `bin/`.

### Example marginalization

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

where the line `1 1 1 ...` is outputted to std::cout whilst the others appear in std::clog (thereby allowing for easy 
redirection of the output).
Each integer on the output line corresponds to the index of the block of vertex v_0, v_1,..., v_n.

Replace `bin/mcmc` with `bin/mcmc_history` to output the state of the system everytime it is sampled.

Option "--use_ppm" enables simpler transition probabilities computation, only possible for probability matrices of the form

    p_in  p_out p_out ...  p_out 
    p_out p_in  p_out ...  p_out 
    p_out p_out p_in  ...  p_out 
      .     .     .    .     .
    p_in  p_out p_out ...  p_in

It must be use in conjunction with `-P p_in p_out` instead of the full matrix

### Example maximization

In the maximization mode, we guess the planted partition by maximizing the likeihood of the partition (with simulated 
annealing).

The call is similar to that of the marginalization mode:

	bin/mcmc -e example_edge_list.txt -P 0.6 0.1 0.1 0.6 -n 20 20 -r -t 1000 --maximize -c exponential	

Both the burn-in and sampling frequency are ignored in the maximization mode.

4 cooling schedules are implemented: `exponential`, `linear`, `logarithmic` and `constant`.

There inverse temperature is given as

    beta(t) = 1/T_0 * alpha^(-t) 	        (Exponential)
    beta(t) = 1/T_0 * [1 - eta * t / T_0]^(-1)  (Linear)
    beta(t) = log(t + d) / c  			(Logarithmic)
    beta(t) = 1 / T_0				(Constant)

where $t$ is the MCMC step. The paramters of these cooling schedule are passed like so:

	-a T_0 alpha    (Exponential)
	-a T_0 eta      (Lienar)
        -a c d          (Logarithmic)
	-a T_0          (Constant)

## Companion article

Please cite

**Finite size analysis of the detectability limit of the stochastic block model**

Jean-Gabriel Young, Patrick Desrosiers, Laurent Hébert-Dufresne, Edward Laurence, Louis J. Dubé

[arXiv:1701.00062](https://arxiv.org/abs/1701.00062)
