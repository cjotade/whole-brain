Fast simulator of the Dynamic Mean Field model of brain dynamics
================================================================

This C++/Matlab software implements the Dynamic Mean Field (DMF) model of brain
dynamics, originally published by Gustavo Deco (UPF) and collaborators (see
references below).

This implementation offers three main advantages over the original Matlab
implementation by Deco *et al*:

1. It is written in C++, which already makes it faster (integration of DMF
   equations is about 3x faster than in Matlab).

2. The Balloon-Windkessel model of BOLD response is run in parallel with the
   DMF model itself, adding a further ~2x speedup when BOLD signals are
   requested.

3. If only BOLD time series are requested, the DMF simulator runs on a circular
   buffer, which radically reduces memory consumption.

Altogether, and to the best of my knowledge, these advantages place this
software amongst the fastest and most efficient fMRI simulators available.

Pedro Mediano, Apr 2020


## Compilation

To compile the code from Matlab, make sure you are in the root folder of
`FastDMF` (so that the Eigen library is in the current folder) and run:

```octave
mex COPTIMFLAGS="-O3" DMF.cpp
```

If for any reason you do not want the DMF and BOLD integrators to run in
parallel (because of system constraints or as a speed test), you can use the
`NO_PARALLEL` compiler flag:

```octave
mex COPTIMFLAGS="-O3" DMF.cpp -DNO_PARALLEL
```


## Usage

The basic usage is as follows:

```octave
bold = DMF(params, nb_steps);
```

Above, `nb_steps` is the number of DMF (not BOLD) steps to integrate, and
`params` is a `struct` with all relevant parameters (see `example.m`). The
resulting matrix, `bold` is of size `NxT`, with N variables and T timesteps.

One can also request firing rates and/or BOLD signals through an optional third
argument to `DMF` (*note*: number of output arguments has to match, otherwise
the code will throw an error):

```octave
rate = DMF(params, nb_steps, 'rate');
bold = DMF(params, nb_steps, 'bold');
[rate, bold] = DMF(params, nb_steps, 'both');
```


## Implementation notes and possible improvements

* The synchronisation between DMF and BOLD integrators is done very crudely by
  spawning and joining BOLD threads every batch. The proper way to do this
  would be with a producer-consumer design (controlled so that DMF and BOLD
  aren't manipulating the same chunk of buffer). Based on the notes below, I
  suspect that for the typical use-case this would have little effect.

* For reasonable values of `batch_size` (>1000) and a 90-region parcellation,
  the BOLD integrator always finishes batches faster than the DMF. This makes a
  funky producer-consumer structure probably unnecessary, although this may not
  hold for smaller models (i.e. with N < 90 ROIs), since DMF scales as O(N^2)
  and BW as O(N).

* The BOLD integrator operates fully component-wise, suggesting that it could
  be further parallelised into array chunks. As mentioned above, however, at
  the moment the simulation is compute-bound by the DMF equations (and not the
  BW), meaning that such a parallelisation won't help much (if at all).

* There are a few points where I think the code might be further optimised for
  a possible extra speedup:

    - Explore the aliasing, in-place operations, and matrix-vector product inside
      the the for-loop in `DMFSimulator::run`.
  
    - Parallelising excitatory and inhibitory computations in the same loop.
  
    - Same considerations for the BW computations in `BOLDIntegrator::compute`.

  Unfortunately, these computations would have pretty short turnaround times,
  so they would ned proper consumer-producer designs (instead of simply
  spawning new threads all the time).


## References

* Deco, G., Hagmann, P., Romani, G. L., Mantini, D. & Corbetta, M. "How local
  excitation-inhibition ratio impacts the whole brain dynamics". J. Neurosci.
  34, 7886–7898, DOI: 10.1523/JNEUROSCI.5068-13.2014 (2014).

* Deco, G., *et al*. "Whole-brain multimodal neuroimaging model using serotonin
  receptor maps explains non-linear functional effects of LSD". Curr. Biol.
  1–10, DOI: 10.1016/j.cub.2018.07.083 (2018).

* Herzog, R., Mediano, P., Rosas, F., Carhart-Harris, R., Sanz Perl, Y.,
  Tagliazucchi, E. & Cofre, R.  "A mechanistic model of the neural entropy
  increase elicited by psychedelic drugs". Pre-print.

