# WIMpy
[**Work in progress**] Python code for doing fits to (mock) dark matter direct detection data.

Details of how to use the code will be added soon!

### Performing fits

The python script `WIM.py` is used to calculate the significance achievable with a particular set of experiments and a particular set of input Dark Matter parameters. It is called as:

`python WIM.py ENSEMBLE MASS INDEX DIR`

Here, `ENSEMBLE` specifies which experimental ensemble to use `A`, `B`, `C` or `D`. `MASS` specifies the input DM mass in GeV. `INDEX` specifies the input DM couplings (in particular, in indexes the points on a grid of couplings).

`WIMpy` generates 50 mock data sets based on the input parameter points. For each data set, it uses a grid-refinement method to calculate the maximum likelihood and therefore the significance for discrimination between Dirac-like and Majorana-like couplings. The significances are output to a file named `Results_pINDEX.txt` in the relative path `DIR`.

### Checking the likelihood calculator

`CompareNgrid.py`

### Plotting
