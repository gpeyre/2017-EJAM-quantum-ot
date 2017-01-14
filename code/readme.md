This toolbox reproduces the numerical illustrations of [our paper](paper/TensorOT.pdf):

G. Peyré, L. Chizat, F-X. Vialard, J. Solomon, [Quantum Optimal Transport for Tensor Field Processing](paper/TensorOT.pdf), Arxiv, 2016

![Example of tensor-valued interpolation](img/interpolation.png)

Functions
-------

The main functions are:
- quantum_sinkhorn: compute an optimal coupling between two tensor fields.
- quantum_interp: use such a coupling to compute an interpolation between the tensor fields.
- quantum_barycenter: compute a barycenter between several fields.

Structure of the toolbox
-------

- test_*.m : the main scripts to reproduce the figures of the article.
- test_simple_*.m: simple tests (display, etc), not related to optimal transport.


Scripts
-------

The main scripts to reproduce the figures of the article are:
- test_2x2_1d.m: 2x2 tensors in 1-D (along a line).
- test_3x3_1d.m: 3x3 tensors in 1-D (along a line).
- test_2x2_2d.m: 2x2 tensors in 2-D (image).
- test_2x2_meshes.m: 2x2 tensors on the tangent plane of a mesh.
- test_2x2_2d_sparse: example of use of the sparse mode (using multiscale computations).

Note that for these algorithms, tensor fields are stored in format (d,d,...) for tensors of size dxd so the first dimensions expose the tensors. This is opposed to the convention used for e.g. vizualization helpers tools.

Helpers
-------

Generic helper functions are in toolbox/.

Helper functions for the quantum-OT are implemented in toolbox_quantum/, in particular:
- load_helpers: load a structure filled with many simple functions.
- LSE: log-sum-exp operator, at the core of Sinkhorn iterations.
- expM / logM: fast matrix exp and log
- plot_tensors_1d/2d: vizualization
- load_tensors_pair: load synthetic 1D and 2D examples.

Helpers functions for mesh processing are in toolbox_geometry/.

Mex
------

You might want to re-compile for your architecture the mex files in tensor_logexp/, which are fast C implementation of logM/expM. The script compile_mex.m should be able to make the compilation.

To enable the use of these fast mex file, set

> global logexp_fast_mode
> logexp_fast_mode = 4; % fast mex

Sparse mode
------

By default, the couplings gamma are explicitly stored as full matrices. To scale to large problem size, it is mandatory to use sparse matrices. This is active if you provide a sparse cost matrix c. This means that you have to "guess" beforehand roughly the location of the transport map. This is usually done by running the code in a coarse to fine way. See test_sparse.m for an example of use.

Thanks
-------

For mesh computations, uses the Matlab toolbox toolbox_connections/ of [Keenan Crane](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/).


Copyright
-------

Copyright (c) 2016 Gabriel Peyré
