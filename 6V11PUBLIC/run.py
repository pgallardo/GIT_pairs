from VGIT import Hypersurfaces
from VGIT import VGIT
from VGIT.OPSubgroups import OPS, OPS_destabilized_pair
from VGIT.Monomials import MonFrozenSet, Monomial
from sympy import Rational
from sympy import *
#Make automatic=0 if you want to be prompted for input.
automatic=1


print "VGIT PROGRAM VERSION",
print VGIT.Version()
print "By Patricio Gallardo and Jesus Martinez-Garcia"
print "VGIT of Hypersurfaces"
#dimension of the variety:
dimension=2
degree=3
debug=0
pairs=1
div_degree=1

if not automatic:
        dimension = input('Enter the dimension of your varieties: ');
        degree = input('Enter the degree of your first variety: ');
        pairs = input('Enter 1 if you want a problem in Variational GIT and 0 if you want a problem in GIT: ');
        if pairs:
                div_degree= input('You have chosen a VGIT problem. Enter the degree of the second variety: ');
        if not (div_degree>=1):
                print 'You chose a non-valid degree. Choosing hyperplane sections as divisors'
                div_degree=1
        debug = input('Debug mode? (Enter 1 for debug mode, 0 for no debug mode: ');

prob=Hypersurfaces.Problem(debug,dimension,degree, pairs, div_degree);
##raw_input('Press Enter to continue')
if debug:
        print('Input received')

#the paramenter gives the t of the GIT problem, which is t=0 at the moment
print 'Solving the problem'
solutions=prob.solve();

if solutions==None:
        print " The solution is empty."
else:
##        Substitute by the following line for latex output.
##        solutions.printout(format_print='latex')
        solutions.printout(pauses=0)

print 'reached the end'
ops=prob.all_1PS()
print ops	
