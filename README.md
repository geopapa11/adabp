# adabp
The Adaptive Belief Propagation (BP) algorithm.

This is research code and builds on the following paper (please cite appropriately):

- Georgios Papachristoudis, John W. Fisher III.
Adaptive Belief Propagation,
ICML 2015.

bibtex: <br>
@inproceedings{ Geopapa2015adabp,<br>
author = {Papachristoudis, Georgios and Fisher III, John W.},<br>
booktitle = {International Conference on Machine Learning (ICML)},<br>
month = {July},<br>
title = {{Adaptive Belief Propagation}},<br>
url = {http://jmlr.org/proceedings/papers/v37/papachristoudis15.pdf},<br>
year = {2015}<br>
}

## Getting Started

Clone this repository:

```
git clone https://github.com/geopapa11/adabp
```

Add the rctreebp folder to the MATLAB path:

```
addpath('rctreebp/');
```

(or you can just go to MATLAB GUI -> Home panel, choose "Set Path" -> "Add Folder" and choose the "rctreebp" folder)

## Main functions

- <tt>AdaBP.m</tt>: Adaptive BP algorithm.
  - Support for discrete forest MRFs.
  - Support for Gaussian forest MRFs.
  - Support for max-product.

## Reproduce paper figures

- If you wish to reproduce Figures 2 and 3 of the ICML paper, please open <tt>adabp_toy.m</tt> and read the help documentation for further instructions.
- If you wish to reproduce Figure 4 of the ICML paper, please open <tt>cpg_toy.m</tt> and read the help documentation for further instructions.
- If you wish to reproduce Figure 5 of the ICML paper, please open <tt>berk_toy.m</tt> and read the help documentation for further instructions.
- If you wish to run a toy example of an application of AdaBP to Gaussian loopy graphs, please type the command <br>
```[mu1, var1, mu2, var2] = adabp_loopy_toy();```

## Comments
It has been noted that RCTreeBP by [Sumer et al., 2011] might crash in Matlab 2014 and above. So, please keep this in mind, when you run the code. If you wish to see the performance of AdaBP alone, you can comment the lines that are related to RCTreeBP in <tt>adabp_toy.m</tt> and <tt>cpg_toy.m</tt>

## Tested on Linux, Mac and Windows

## Author
Giorgos Papachristoudis (geopapa)
