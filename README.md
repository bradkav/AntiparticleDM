# AntiparticleDM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.815457.svg)](https://doi.org/10.5281/zenodo.815457)

Python code for calculating the prospects of future direct detection experiments to discrimination between Majorana and Dirac Dark Matter (i.e. to determine whether Dark Matter is its own antiparticle).

Direct detection event rates and mock data generation are taken care of by a variation of the `WIMpy` code - available [here](https://github.com/bradkav/WIMpy/tree/Antiparticle).

With this code, the results of arXiv:XXXX.XXXXX should be **entirely reproducible**. Follow the instructions under [below](#repro) if you want to reproduce those results. If you find any mistakes or have any trouble reproducing the results, please get in touch.

If you have any questions, comments, bug-reports etc., please contact Bradley Kavanagh at bradkav@gmail.com. 

### Version History

**Version 1.0 (23/06/2017):** Initial release, including all results and plots from the paper.

## Contents

- `calc`: core code for calculating the statistical significance for discriminating between Dirac and Majorana Dark Matter (DM).
- `scripts`: scripts for reproducing results from the paper (NB: some may need to be implemented on a computing cluster...)
- `analysis`: scripts for processing the results and generating plots.
- `results`: data products for a range of DM masses, couplings and experimental ensembles.
- `plots`: plots from arXiv:XXXX.XXXXX (and others).

## Reproducing the results <a name="repro"></a>

The majority of the code is written in `python`, and requires the standard `numpy` and `scipy` libraries. For plotting, `matplotlib` is also required.

### Performing likelihood fits

***NEED TO MENTION THE BINDER***

Code for generating mock data sets and performing likelihood fits are found in the `calc` folder. The module `CalcDiscrimination.py` does the heavy lifting. If you want to calculate the discrimination significance for a particular ensemble and parameter point, use:

```python
import CalcDiscrimination as CD
CD.CalcDiscrim(ensemble, m0, f, r_np, outfile=None)
```

where `ensemble` specifies which experimental ensemble you consider (A, B, C or D), `m0` is the DM mass in GeV, `f` is the interference factor (Eq. 6) and `r_np` is the ratio of neutron to proton couplings (Eq. 5 over Eq. 4). You can also specify an file to save the results to (`outfile`). 

**Note** that the number of mock data samples (`N_samples`) to take is specified in the file `params.txt`. For illustration purposes, this is set to `N_samples = 10`. However, we recommend using a value of `N_samples = 100` to get enough statistics.

For calculating the discrimination significance over a grid of the input couplings, you can run `RunFits_couplings.sh`. 

For calculating the discrimination significance as a function of exposure (for a fixed input), you can `RunFits_exposure.sh`.

Note that these scripts will take a long time to run (think hours to days...). In practice then, you'll probably want to run things on a computing cluster. For this, we provide two python files `RunMPI_couplings.py` and `RunMPI_exposure.py`, which are MPI-enabled and take care of running large numbers of fits in parallel. To use these, `mpi4py` is required.
<!---
For each data set, it uses a grid-refinement method to calculate the maximum likelihood and therefore the significance for discrimination between Dirac-like and Majorana-like couplings. The significances are output to a file named `Results_pINDEX.txt` in the relative path `DIR`.
--->

### Generating plots

Scripts for generating plots from the results are in the `analysis/` folder. To (re-)generate all the plots from the paper, simply run `analysis/GeneratePlots.sh`.

### Checking the likelihood calculator

`CompareNgrid.py`





## Citation

If you make use of the code or the numerical results, please cite the project as:

> Kavanagh, B. J., Queiroz, F. S., Rodejohann, W., Yaguna, C. E., "AntiparticleDM" (2017), https://github.com/bradkav/AntiparticleDM/, [doi:10.5281/zenodo.815457](http://dx.doi.org/10.5281/zenodo.815457)

Please also cite the associated papers:

> Kavanagh, B. J., Queiroz, F. S., Rodejohann, W., Yaguna, C. E., "Prospects for determining the particle/antiparticle nature of WIMP dark matter with direct detection experiments" (2017), arXiv:XXXX.XXXXX

> Queiroz, F. S., Rodejohann, W., Yaguna, C. E., "Is the dark matter particle its own antiparticle?", Phys. Rev. D 95 (2017) 095010, [arXiv:1610.06581](https://arxiv.org/abs/arXiv:1610.06581)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
