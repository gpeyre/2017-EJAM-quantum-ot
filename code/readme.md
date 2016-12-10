This code implements quantum optimal transport (transport and barycenters) using an entropic regularization (corresponding to matrix-valued generalization of Sinkhorn's algorithm).

The core algorithms are:
- quantum_barycenters
- quantum_sinkhorn
Note that for these algorithm, tensor fields are stored in format (2,2,...) so the first dimensions expose the tensors. This is opposed to the convention used for e.g. vizualization helpers tools.

Helper functions are implemented in toolbox_quantum/, in particular:
- load_helpers: load a structure filled with many simple functions.
- LSE: log-sum-exp operator, at the core of Sinkhorn iterations.
- expM / logM: fast matrix exp and log
- plot_tensors_1d/2d: vizualization
- load_tensors_pair: load synthetic 1D and 2D examples.

Copyright (c) 2016 Gabriel Peyre
