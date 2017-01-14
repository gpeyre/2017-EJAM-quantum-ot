-----------------------------------------------------------
Example Code for "Trivial Connections on Discrete Surfaces"
-----------------------------------------------------------
Keenan Crane
June 15, 2010

This program computes smooth direction fields on meshes using the method
described in Crane et al, "Trivial Connections on Discrete Surfaces."
The primary purpose of this code is to illustrate how the algorithm works;
it is likely not suitable for large meshes.  Furthermore, for simplicity,
this version of the code works only on genus-0 surfaces without boundary,
i.e., topological spheres.

You can test the code from the MATLAB command line by simply typing 'run'
from the install directory.  You should be presented with a plot that looks
like the one in test.png, exhibiting a single singularity of index +2 on
the sphere.

This code was developed for an introductory undergraduate course on discrete
differential geometry, and follows three prior assignments:

  Topological Invariants of Discrete Surfaces
  http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw1.pdf

  Mesh Smoothing
  http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw2.pdf

  Vector Field Decomposition
  http://www.cs.caltech.edu/~keenan/ddg_exercises/ddg_hw2.pdf

The material above builds all of the necessary concepts (both code and math)
from the ground up, and should provide a concrete understanding of most of
the ideas used in this code -- the remaining ideas are discussed in the paper

   "Trivial Connections on Discrete Surfaces"
   Keenan Crane, Mathieu Desbrun, Peter Schröder
   SGP 2010 / Computer Graphics Forum
   http://www.cs.caltech.edu/~keenan/pdf/connections.pdf

Finally, please do not hesitate to send bugs or questions to keenan@cs.caltech.edu

---------------------
 LICENSE
---------------------

The following terms and conditions apply to the materials contained in this
archive -- they are based on the FreeBSD open source software license:

Copyright 2010 Keenan Crane. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the author and should not be interpreted as representing official policies,
either expressed or implied, of any other person or institution.

