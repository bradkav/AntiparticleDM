## AntiparticleDM - Calc

#### Running the core of the code

The module `CalcDiscrimination.py` does the heavy lifting. If you want to calculate the discrimination significance for a particular ensemble and parameter point, use:

```python
import CalcDiscrimination as CD
CD.CalcDiscrim(ensemble, m0, f, r_np, outfile=None)
```

where `ensemble` specifies which experimental ensemble you consider (A, B, C or D), `m0` is the DM mass in GeV, `f` is the interference factor (Eq. 6) and `r_np` is the ratio of neutron to proton couplings (Eq. 5 over Eq. 4). You can also specify a file to save the results to (`outfile`). 

**Note:** the number of mock data samples (`N_samples`) to take is specified in the file `params.txt`. For illustration purposes, this is set to `N_samples = 10`. However, we recommend using a value of `N_samples = 100` to get enough statistics.

#### How it works

The files in `DDexpt` specify the properties of a number of direct detection experiments (numbers of protons, neutrons, exposure, mass fractions, etc.), which are loaded by the `WIMpy` code. The `WIMpy` code also calculates event rates, generates mock data and pre-tabulates some of the likelihood information.

The `CalcLikelihood.py` then takes care of finding the maximum likelihood under the Majorana and Dirac DM hypotheses, using a grid-refinement method. In practice, this just means rescaling the number of expected events observed at each value of the couplings and calculating the likelihood from this.

For a given sample, the code then calculates the significance with which a Majorana DM particle can be excluded. This procedure is repeated for a large number of samples.

#### Useful files

- `CalcDisc-vs-Couplings.py`: Calculate discrimination significance as a function of 'coupling index' (which specifies the position on a grid of the DM-nucleon couplings.
- `CalcDisc-vs-Exposure.py`: Calculate discrimination significance as a function of 'exposure index' (which specifies the exposure of the experiments other than Xe and Ar). 
- `CompareNgrid.py`: Compare the likelihood calculated with different numbers of grid points to show good convergence.

Check the files themselves to see how they work in a bit more detail. You can also check the `scripts` folder to see how you would use them in practice.
