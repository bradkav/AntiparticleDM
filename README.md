# AntiparticleDM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.815457.svg)](https://doi.org/10.5281/zenodo.815457)

Python code for calculating the prospects of future direct detection experiments to discrimination between Majorana and Dirac Dark Matter (i.e. to determine whether Dark Matter is its own antiparticle).

Direct detection event rates and mock data generation are taken care of by a variation of the `WIMpy` code - available [here](https://github.com/bradkav/WIMpy/tree/Antiparticle).

With this code, the results of arXiv:XXXX.XXXXX should be **entirely reproducible**. Follow the instructions under [below](#repro) if you want to reproduce those results. If you find any mistakes or have any trouble reproducing the results, please get in touch.

If you have any questions, comments, bug-reports etc., please contact Bradley Kavanagh at bradkav@gmail.com. 

### Version History

## Contents

## Reproducing the results <a name="repro"></a>

### Getting started

### Performing fits

The python script `CalcDiscrimination.py` is used to calculate the significance achievable with a particular set of experiments and a particular set of input Dark Matter parameters. It is called as:

`python CalcDiscrimination.py ENSEMBLE MASS INDEX DIR`

Here, `ENSEMBLE` specifies which experimental ensemble to use `A`, `B`, `C` or `D`. `MASS` specifies the input DM mass in GeV. `INDEX` specifies the input DM couplings (in particular, in indexes the points on a grid of couplings).

`WIMpy` generates 50 mock data sets based on the input parameter points. For each data set, it uses a grid-refinement method to calculate the maximum likelihood and therefore the significance for discrimination between Dirac-like and Majorana-like couplings. The significances are output to a file named `Results_pINDEX.txt` in the relative path `DIR`.

### Checking the likelihood calculator

`CompareNgrid.py`

### Plotting



## Citation

If you make use of the code or the numerical results, please cite the project as:

`Kavanagh, B. J., Queiroz, F. S., Rodejohann, W., Yaguna, C. E., AntiparticleDM (2017), https://github.com/bradkav/AntiparticleDM/`

Please also cite...

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
