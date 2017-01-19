# GIT_pairs

Varitions of GIT quotients package

This page is a complement to the articles Moduli of cubic surfaces and their anticanonical divisors and Variations of geometric invariant quotients for pairs, a computational approach authored by Patricio Gallardo and Jesus Martinez-Garcia
(see home page). The following instructions are intended for the readers that want to use the raw output of our program, and for the readers that want to run the program by themselves.

License

The package VGIT-0.6.12, was created to complement these articles and we expect it can be used by other mathematicians to solve problems. The code is based on the theory developed in the article Variations of geometric invariant quotients for pairs, a computational approach and a rough idea of the algorithm can be found there. More detailed algorithms will appear in an upcoming article, together with further applications. The source code and data, but not the text of this article, are released under a Creative Commons CC BY-SA 4.0 license. If you make use of the source code and/or data in an academic or commercial context, you should acknowledge this by including a reference or citation to Variations of geometric invariant quotients for pairs, a computational approach —in the case of the code— or to Moduli of cubic surfaces and their anticanonical divisors —in the case of the data for cubic surfaces. 

Description of the output.

We describe the file so that it can be understood in the context our article. For a description of the algorithms that generate the output, we refer the reader to our article Variation of geometric invariant quotients for pairs a computational approach.

1. The program returns a list $t_i$ including all GIT walls and a representative for each GIT chamber.
2. The program returns the fundamental set of one-parameter subgroups. The program also tells you which one-parameter subgroups are not necessary for the GIT analysis.
3. For each $t_i$, the program returns a list of one-parameter subgroups that generate non-stable maximal sets. The program recognizes the changes in the non-stable maximal sets with respect the previous $t_{i-1}$. There are three types of outputs:

 3.1. If a non-stable set of monomial has appeared for an earlier $t_i$, then it is denoted as "same".
 
 3.2. If it becomes semistable, then it is denoted as "same but becomes semistable."
 
 3.3. If it is a new set of monomials, then it is denoted as "new".

4. Finally, the program tells the reader if a maximal non-stable set has disappeared when to compare it with the previous $t_{i-1}$. For each one-parameter subgroup that it is not denoted as "same", the program returns the monomials in the maximal non-stable sets.
5. The program recognizes if the configuration is semi-stable or not. In the former case, it returns the potential close orbit of that pairs.

How to run the code.

We recommend the reader to sign up and use https://cloud.sagemath.com. Otherwise, it will be necessary to install Python 2.7 and the required libraries such as "Sympy". The raw code can be downloaded here or at http://guests.mpim-bonn.mpg.de/martinezgarcia/code.html

Description of the package

For a precise description of the coded algorithms, we refer the reader to our article "Variation of geometric invariant quotients for pairs a computational approach." The code in the folder VGIT has the following files:

1. The file " __init__.py" is required to make Python treat the directories as containing packages. Then, we can import the modules 'VGIT', 'Hypersurfaces', 'Monomials', and 'OPSubgroups'.

2. The file "Monomials.py" contains the sub-package Monomials with the class Monomials to deal with one monomial and the class MonFrozenSet to deal with a set of monomials. In this file, the reader can find the Hilbert_Mumford function and the function that generates the normalized one-parameter subgroups from a list of monomials.

3. The file OPSubgroups.py contains the sub-package OPSubgroups with the class OPS and all the methods to work with them. It imports from the class "Monomials", "associated1ps" and "MonFrozenSet"

4. The file "VGIT.py" contains the metaclasses "Problem", "Solution" and "Solution_t" that defines the VGIT problem and it stores the solutions.

5. The file "Hypersurface.py" contains the classes "Problem" and VGIT.Hypersurfaces.Solution_t". We store the VGIT problem in the class "Problem." We store the solution to our VGIT problem in the class "VGIT.Hypersurfaces.Solution_t."

6. In the file run.py, the reader can change the dimension, and degree of the VGIT problem.
