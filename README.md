# GIT_pairs

This package was created to complement the article  "VGIT for pairs, a computational approach" which is joint work with Jesus Martinez-Garcia  (see http://arxiv.org/abs/1602.05282). We expect our software can be used by other mathematicians to solve problems related to VGIT. The source code and data are released under a Creative Commons CC BY-SA 4.0 license. If you make use of the code and data in an academic or commercial context, you should acknowledge a reference or citation to our work.   

The code in the folder VGIT has the following structure:

1. The file " __ init __.py" is required to make Python treat the directories as containing packages. Then, we can import the modules 'VGIT', 'Hypersurfaces', 'Monomials', and 'OPSubgroups'.
2. The file "Monomials.py" contains the sub-package Monomials with the class Monomials to deal with one monomial and the class MonFrozenSet to deal with a set of monomials. In this file, the reader can find the Hilbert_Mumford function and the function that generates the normalized one-parameter subgroups from a list of monomials.
3. The file OPSubgroups.py contains the sub-package OPSubgroups with the class OPS and all the methods to work with them. It imports from the class "Monomials", "associated1ps" and "MonFrozenSet"
4. The file "VGIT.py" contains the metaclasses "Problem", "Solution" and "Solution_t" that defines the VGIT problem and it stores the solutions.
5. The file "Hypersurface.py" contains the classes "Problem" and VGIT.Hypersurfaces.Solution_t.  We store the VGIT problem in the class "Problem." We store the solution to our VGIT problem in the class "VGIT.Hypersurfaces.Solution_t."
6. In the file run.py, the reader can change the dimension, and degree of the VGIT problem.
