Adaptive inference describes the task of recomputing the results of an
inference problem which is slightly different than one which has been already
solved. It uses concepts from incremental computing, in which the amount of
work we expect to do to solve a problem is (hopefully) proportional to the
amount of change in the problem from our already-solved version. In particular,
we use rake and compress trees (RC-Trees) to cache sufficient statistics of the
inference process in such a way that a large fraction of them are guaranteed to
be reusable no matter what (small) changes are made to the original graphical
model, including changing factors or potential functions, adding or removing
variables, or modifying the graph connectivity. For a system of N variables,
the sum-product version of the code can update the marginal distribution of K
specified variables in O(K log N/K) time. The max-product version of the code
keeps track of an optimal configuration of variables (the MPE or MAP
configuration), and after any change which induces K changes in the optimal
configuration, finds all those changes in O(K log N/K) time, without the
locations needing to be specified.

For more information, please see our publications on the subject:

    * Adaptive Bayesian Inference, NIPS 2007
    * Adaptive Inference in General Graphical Models, UAI 2008
    * Adaptive Updates for MAP Configurations with Applications to Bioinformatics, SSP 2009 


==============
Using the Code
==============

This code is mainly intended for educational purposes, to provide a working
implementation for anyone wishing to make use of the algorithms in their own
work. It is not optimized, nor does it have the greatest of user interfaces.

The Matlab functions given here were used in the NIPS 2007 publication, and
assume that the underlying graph is tree-structured. The SSP 2009 paper also
used mainly tree-structured graphs, including HMMs, and the code should work
for these problems as well, but I haven't tried it.

The code includes versions for sum-product (computing marginal distributions)
and dynamic programming for optimization (max-product to compute sufficient
statistics, then selecting a configuration of the variables in the downward
pass). The functions corresponding to each are designated "SP" and "MP",
respectively. 

=======================
Using the MEX Functions
=======================

The algorithms involved in this work rely heavily on pointers and memory
re-allocation, neither of which is a strong point of the Matlab programming
environment. To make the algorithms scale in the proper way, certain aspects of
Matlab's memory management system were bypassed; this is done using MEX files
to directly access and change the array and cell array contents. However, these
functions are undocumented and are not guaranteed; they worked for me but may
not work for you. Here is the disclaimer in the C files:

WARNING! WARNING! WARNING! This MEX file is not for the faint of heart. It
attempts to circumvent Matlab's memory-handling functions in order to give more
efficient code, and uses undocumented functions to do so. It comes with no
guarantees whatsoever. It may not work on other versions of Matlab, or other
platforms, or even at all. It may result in memory leaks, segmentation faults,
destroy data, set your computer on fire, or suck it into a black hole, for all
I know. If you do not accept these possible risks, do not use this MEX file.

If that doesn't deter you, assuming you have MEX properly set up (see my KDE
toolbox page for some pointers on that subject) you should be able to compile
the functions using the "buildmex" script. 

====================
 COPYRIGHT / LICENSE
====================

The adaptive inference code for Matlab was written by Alexander Ihler, and are
copyrighted under the (lesser) GPL:

    Copyright (C) 2007-2008 Alexander Ihler

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; version 2.1 or later. This program is distributed in the
hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details. You should have received a copy
of the GNU Lesser General Public License along with this program; if not, write
to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

The authors may be contacted via email at: ihler (at) alum (.) mit (.) edu 

=======
Updates
=======

See the authors' webpages for the most current version of the code and documentation.
