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



