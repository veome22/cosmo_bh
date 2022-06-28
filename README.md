# cosmo_bh
Warning: This is a work in progress!

In this project, we simulate a set of Stellar Binaries in isolation using [COMPAS](https://arxiv.org/abs/2109.10352), with the intention of studying the connection of Compact Object (CO) binary mergers to their host galaxy properties.

The range of initial Stellar Binaries simulated is [5e1, 5e2, 5e3, 5e4, 5e5, 5e6].

The range of metallicities is [0.0001, 0.00018847, 0.0003552, 0.00066943, 0.00126166, 0.00237782, 0.0044814, 0.00844598, 0.01591789, 0.03] $Z_{abs}$.

The notebook containing the results of the simulations is [bh_properties.ipynb](../main/bh_properties.ipynb).


# Drawing a Sample of BBHs at a given Star Forming Mass, Metallicity
This repo uses the tools defined in the [cosmo_bbh_tools](../main/cosmo_bbh_tools/) directory to sample a population of BBHs using any Star Forming Mass and Metallicity value, as long as it lies within the simulated range. The main function to sample a population is called [sample_bbh_from_sfm_met()](../main/cosmo_bbh_tools/cosmo_bbh_tools.py#L298-L347).

A complete example of how to sample a BBH population, along with code to load the relavant data, is shown in [draw_bbhs.ipynb](../main/draw_bbhs.ipynb). 
**PLEASE START HERE.**


