##################################################################
#Package name: VGIT.VGIT
#Creation date: 14/05/20
#Description:
#   Contains the metaclasses Problem and Solution that defines VGIT
#       in several contexts, as well as several auxiliary functions
##################################################################

########################
#GLOBAL VARIABLES
########################

version="0.6.11"

########################
#AUXILIARY FUNCTIONS
########################

def Version():
    return version

########################
#CLASES
########################

########################
#Class: Problem
#Attributes:
#       debug: 1 for debug mode
#       dim: dimension of the objects in the problem
#           e.g. a problem of surfaces, dim=2
#Description: Metaclass to define a VGIT problem
#Methods:
#   __init__(self,debug,dimension): creates new objects
#   solve_t(self,t): solves the Problem for fixed rational t
########################
class Problem(object): 
	def __init__(self,debug,dimension):
            self.dim=dimension
            self.debug=debug
            if debug==1:
                print 'debug ',self.debug, 'dimension ', self.dim
        def solve(self):
            return self.solve_t(0)

        def solve_t(self,t):
            sol=Solution_t(self.debug, self.dim);
            if self.debug==1:
                print 'Method Problem.solve_t not yet implemented.'
            return sol
            
########################
#Class: Solutions
#Attributes:
#       debug: 1 for debug mode
#       dimension: dimension of the objects in the problem
#           e.g. a problem of surfaces, dim=2
#Description: Metaclass to store the solution to a VGIT problem
#Methods:
#   __init__(self, debug,dimension): creates a new solution
#   print_sol(self): prints the information of the solution
########################
class Solutions(object):
        def __init__(self, debug, dimension):
            self.debug=debug
            self.dimension=dimension

            
            
########################
#Class: Solution_t
#Attributes:
#       debug: 1 for debug mode
#       dimension: dimension of the objects in the problem
#           e.g. a problem of surfaces, dim=2
#Description: Metaclass to store the solution to a VGIT problem
#Methods:
#   __init__(self, debug,dimension): creates a new solution
#   print_sol(self): prints the information of the solution
########################

class Solution_t(object): # VGIT problem
	def __init__(self):
            pass

        def print_sol(self):
            print 'Dimension: ', self.dim
            if self.debug==1:
                    print 'Method Solution_t.print not yet implemented.'
            
