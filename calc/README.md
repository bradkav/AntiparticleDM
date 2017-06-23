## AntiparticleDM - Calc

The module `CalcDiscrimination.py` does the heavy lifting. If you want to calculate the discrimination significance for a particular ensemble and parameter point, use:

```python
import CalcDiscrimination as CD
CD.CalcDiscrim(ensemble, m0, f, r_np, outfile=None)
```

where `ensemble` specifies which experimental ensemble you consider (A, B, C or D), `m0` is the DM mass in GeV, `f` is the interference factor (Eq. 6) and `r_np` is the ratio of neutron to proton couplings (Eq. 5 over Eq. 4). You can also specify an file to save the results to (`outfile`). 

**Note** that the number of mock data samples (`N_samples`) to take is specified in the file `params.txt`. For illustration purposes, this is set to `N_samples = 10`. However, we recommend using a value of `N_samples = 100` to get enough statistics.
