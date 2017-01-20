##################################################################
#Package name: VGIT.Hypersurfaces
#Creation date: 14/05/20
#Description: Contains the subpackage Hypersurfaces, for VGIT 
#   problems about hypersurfaces
##################################################################
from VGIT import Problem, Solution_t, Solutions
from Monomials import MonFrozenSet, construct_boundary
from OPSubgroups import OPS_set, OPS, OPS_destabilized_pair, semistable
import time
from sympy import Rational
########################
#Class: VGIT.Hypersurfaces.Problem
#Parent class: VGIT.Problem
#Attributes:
#       debug: 1 for debug mode
#       dimension: dimension of the problem, e.g. for surfaces dimension=2
#       deg: degree of the problem.
#Description: Class to define a VGIT problem of hypersurfaces
#Methods:
#   __init__(self,debug,dimension): creates new objects
########################
def compute_t(dim, deg, all1PS, allMonomials, allDivisors):
        final=Rational(deg, dim+1)
        t_set=set([Rational(0,1), final])
        if (all1PS==None or allMonomials==None or allDivisors==None):
                return None
        for gamma in all1PS:
                for variety in allMonomials:
                        for divisor in allDivisors:
                                if (gamma.Prod(divisor)!=0):
                                        t=Rational(-gamma.Prod(variety), gamma.Prod(divisor))
                                else:
                                        break
                                if (t<final and t>0):     #Eliminates those we already know have all pairs unstable.
                                        t_set.add(t)
        t_list=sorted(list(t_set))
        return t_list

class Problem(Problem):
	def __init__(self,debug,dimension,degree, pairs=0, div_degree=1):
		super(Problem, self).__init__(debug,dimension)
		self.deg=degree;
		self.pairs=pairs;
		self.div_degree=div_degree;
		self.all1PS=None;
		#Since the monomials are in the ambient projective space
		#we must add 2 to the dimension (one to make it projective,
		#an andother for the fact that X is hypersurface of P^{n+1}
		if debug: print 'Creating set with all monomials.'
		self.allMonomials=MonFrozenSet(self.dim+2, self.deg, 0)
		self.boundaryOPS=MonFrozenSet(self.dim+2,2,construct_boundary(self.dim+2))
		if pairs:
                        self.allPairs=MonFrozenSet(self.dim+2,self.deg,[], gamma=None,pairs=self.pairs, div_degree=self.div_degree,t=0)
			self.allDivisors=MonFrozenSet(self.dim+2, self.div_degree, 0)
                self.all1PS=OPS_set(self.dim+2, self.allMonomials, self.boundaryOPS, self.pairs, self.allDivisors)
		if pairs:
			self.t_list=compute_t(self.dim, self.deg, self.all1PS, self.allMonomials, self.allDivisors)
		else:
                        self.allPairs=None
                        self.allDivisors=None

        def t_s(self):
                return self.t_list
	def all_pairs(self):
                if (not self.pairs):
                        return None
                else:
                        return self.allPairs
        def all_monomials(self):
                return self.allMonomials
        def all_divisors(self):
                return self.allDivisors
        def all_1PS(self):
                if (self.all1PS):
                        return self.all1PS
                else:
                        return None                                        
	def solve(self):
                solutions=Solutions(self.debug,self.dim, self.deg, self.pairs, self.allPairs,self.allMonomials,self.allDivisors)
                for i in range(len(self.t_list)):
                        t=self.t_list[i]
                        if i!=0:
                                solution_t0=solution_t
                        solution_t=self.solve_t(t)
                        if (solution_t):
                                solution_t.type_wall()
                                if i!=0:
                                        if solution_t0!=solution_t:
                                                solutions.add_solution(t, solution_t)
                                        else:
                                                solutions.add_discarded_solution(t)
                                else:
                                        solutions.add_solution(t, solution_t)
                        else:
                                return None
                        if i!=len(self.t_list)-1:
                                solution_t0=solution_t
                                t=(self.t_list[i]+self.t_list[i+1])/2
                                solution_t=self.solve_t(t)
                                if solution_t:
                                        solution_t.type_chamber()
                                        if solution_t0!=solution_t:
                                                solutions.add_solution(t, solution_t)
                                        else:
                                                solutions.add_discarded_solution(t)
                                else:
                                    return None
                return solutions  

                        
	def solve_t(self,t):
		sol=Solution_t(t, self.allPairs, self.allMonomials, self.allDivisors)
		if self.debug:
			print 'Creating set with all 1PS.'
		if self.debug:
			print 'It found all 1PS'                

		all_M_sets=set()
		if self.pairs:
			if self.debug: print 'Solving problem for pairs.\n Creating a set of all sets of pairs destabilised by a 1-parameter subgroup'
			dict_M_gamma={}
			for gamma in self.all1PS:
                                dict_M_gamma[gamma]=[]
                                for boundary in self.allDivisors:
                                        M_gamma_pair=OPS_destabilized_pair(self.dim+2, self.deg, gamma, t, self.allMonomials, self.allDivisors, boundary)
                                        dict_M_gamma[gamma].append(M_gamma_pair)
                                        all_M_sets.add(M_gamma_pair)
			if self.debug: print 'Finding maximal sets of monomials N^+(lambda)'
			for gamma in self.all1PS:
				M_gamma_list=dict_M_gamma[gamma]
				for M_gamma_pair in M_gamma_list:
                                        flag=1
                                        Mgamma_variety, Mgamma_divisor = M_gamma_pair
                                        for M in all_M_sets:
                                                M_variety, M_divisor = M
                                                if (((Mgamma_variety<=M_variety) and (Mgamma_divisor<M_divisor)) or ((Mgamma_variety<M_variety) and (Mgamma_divisor<=M_divisor))) :
                                                        flag=0
                                                        break
                                        if flag:
                                                if semistable(self.dim+2, self.deg, t, M_gamma_pair):
                                                        sol.add_non_stable(gamma,M_gamma_pair, 1)
                                                        sol.add_closed_orbit(gamma, gamma.annihilator(self.dim+2, self.deg,M_gamma_pair, t))
                                                else:
                                                        sol.add_non_stable(gamma,M_gamma_pair, 0)
                                                
			if self.debug: print 'Problem solved'
		else:   
			if self.debug: print 'Solving problem for varieties.\n Creating a set of all sets of monomials destabilised by a 1-parameter subgroup'
			for gamma in self.all1PS:
				M_gamma=MonFrozenSet(self.dim+2,self.deg,[],gamma,self.pairs, self.div_degree,t)
				all_M_sets.add(M_gamma)
			if self.debug: print 'Finding maximal sets of monomials N^+(lambda)'
			for gamma in self.all1PS:
				M_gamma=MonFrozenSet(self.dim+2, self.deg,[], gamma,self.pairs, self.div_degree,t)
				flag=1
				for M in all_M_sets:
					if M_gamma<M:
						flag=0
						break
				if flag:
					sol.add_non_stable(gamma,M_gamma) #***ADD: SS OR NOT SS
			if self.debug: print 'Problem solved'

			
			
	       
		return sol
	def allPairs_list(self):
	    lista= list(self.allPairs)
	    return lista


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
class Solutions(Solutions):
       	def __init__(self,debug,dimension, degree, pairs=0, allPairs=None,allMonomials=None,allDivisors=None):
                super(Solutions, self).__init__(debug,dimension)
                self.deg=degree
                self.pairs=pairs
                self.t_list=[]
                self.walls=[]
                self.chambers=[]
                self.discarded_t_list=[]
                self.solutions={}
                self.allPairs=allPairs
                self.allMonomials=allMonomials
                self.allDivisors=allDivisors
        def add_discarded_solution(self, t):
                self.discarded_t_list.append(t)
        def add_solution(self, t,solution_t):
                self.t_list.append(t)
                self.solutions[t]=solution_t
                if solution_t.t_type()=='chamber':
                        self.chambers.append(t)
                else:
                        self.walls.append(t)
                return
        def printout(self, format_print=None, pauses=1):
		if self.pairs:
			print('\n\nThis is a problem of VGIT, we parametrize pairs of hypersurfaces in projective space and divisors given by restriction of hypersurfaces in that projective space')
		else:
			print('\n\nThis is a problem of GIT, we parametrise hypersurfaces in projective space')
		print 'Dimension: ', self.dimension
		print 'Degree: ', self.deg
		if self.pairs:
			print 'The are %d walls, including the first and last' % len(self.walls)
			print 'There are %d chambers' % len(self.chambers)
			print 'The walls are:'
			print self.walls
			print 'The chambers are:'
			print self.chambers
			print 'Both walls and chambers:'
			print self.t_list
##			print 'The following are \'false positives\' that have been discarded'
##			print self.discarded_t_list
			if pauses:
                                raw_input('Press enter to continue')
			print '\n\n'
		for t in self.t_list:
                        self.solutions[t].printout(self.pairs, self.allPairs, format_print, pauses)
                
	

########################
#Class: VGIT.Hypersurfaces.Solution_t
#Parent class: VGIT.Solution_t
#Attributes:
#       debug: 1 for debug mode
#       dimension: dimension of the problem
#Description: Class to store the solution to a VGIT problem
#               of hypersurfaces
#Methods:
#   __init__(self, debug,dimension): creates a new solution
#######################
class Solution_t(Solution_t): 
	def __init__(self,t=0, allPairs=None, allMonomials=None, allDivisors=None, sol_type='wall'):
		super(Solution_t, self).__init__()
		self.t=t
		self.non_stable={}
		self.non_stable_set=set()
		self.semistable={}
		self.semistable_set=set()
		self.closed_orbit={}
		self.closed_orbit_set=set()
		self.type=sol_type
		return None
	def __eq__(self, other):
                if (len(self.non_stable_set)!=len(other.non_stable_set)) or (len(self.semistable_set)!=len(other.semistable_set)):
                        return False
                for non_stable_pair in self.non_stable_set:
                        if non_stable_pair not in other.non_stable_set:
                                return False
                for semistable_pair in self.semistable_set:
                        if semistable_pair not in other.semistable_set:
                                return False
                return True
	def __ne__(self, other):
                return (not self==other)

        def add_non_stable(self, gamma, M_gamma, semistable=1):
		#First we check that the set M_gamma is not in the solutions as
		#N^+ of a different 1-parameter subgroup
                if M_gamma in self.non_stable_set:
                        return
                gamma_list=list()
                if (gamma in self.non_stable.keys()):
                        gamma_list=self.non_stable[gamma]+[M_gamma]
                else:
                        gamma_list.append(M_gamma)
                self.non_stable[gamma]=gamma_list
		self.non_stable_set.add(M_gamma)
		if semistable:
                        self.add_semistable(gamma, M_gamma)
		return
        def add_closed_orbit(self, gamma, M_gamma):
		#First we check that the set M_gamma is not in the solutions as
		#N^+ of a different 1-parameter subgroup
                if M_gamma in self.closed_orbit:
                        return
                gamma_list=list()
                if (gamma in self.closed_orbit.keys()):
                        gamma_list=self.closed_orbit[gamma]+[M_gamma]
                else:
                        gamma_list.append(M_gamma)
                self.closed_orbit[gamma]=gamma_list
		self.closed_orbit_set.add(M_gamma)
		return

	def add_semistable(self, gamma, M_gamma):
		#First we check that the set M_gamma is not in the solutions as
		#N^+ of a different 1-parameter subgroup
                if M_gamma in self.semistable_set:
                        return
                gamma_list=list()
                if (gamma in self.semistable.keys()):
                        gamma_list=self.semistable[gamma]+[M_gamma]
                else:
                        gamma_list.append(M_gamma)
                self.semistable[gamma]=gamma_list
		self.semistable_set.add(M_gamma)
		return
	def OPS_list(self):
		return list(self.non_stable.keys())
	def OPS_list_semistable(self):
		return list(self.semistable.keys())
	def pairs_non_stable_set(self):
                return self.non_stable_set
	def pairs_semistable_set(self):
                return self.semistable_set
	def pairs_non_stable_list(self):
                return list(self.non_stable_set)
        def pairs_semistable_list(self):
                return list(self.semistable_set)
        def type_wall(self):
                self.type='wall'
                return
        def type_chamber(self):
                self.type='chamber'
        def t_type(self):
                return self.type
        def is_semistable(self, M_gamma):
                if M_gamma in self.semistable_set:
                        return 1
                else:
                        return 0
                


	def printout(self, pairs, allPairs, format_print=None, pauses=1):
		if pairs:
			print 'Solution for t=', self.t,
			print ' which is a %s.' % self.type

		print 'There are %d non-stable maximal sets of monomials' % len(self.non_stable_set),
		print 'of which, %d are semistable' % len(self.semistable_set)
		print 'We have selected %d destabilising test configurations' % len(self.non_stable.keys())
		destabilising_OPS=list(self.non_stable.keys())
		print 'List of destabilising 1-parameter subgroups'
		for gamma in destabilising_OPS:
			print gamma
		if pauses:
                        raw_input('Press enter to continue')
		print '\nList of 1-parameter subgroups and their corresponding set of monomials which are not stable:'
		print '+++++++++++++++++++++++++++++++++++++++++++++++\n\n'
		i=1
		if pairs:
			for gamma, M_gammas_list in self.non_stable.iteritems():
                                print '(',i,') ', 'OPS=',gamma,'\n\n'
                                for M_gamma in M_gammas_list:
                                        variety, divisor=M_gamma
                                        print 'N^+(',gamma,')'
                                        print 'Monomials variety (general element):',
                                        for monomial in variety:
                                                if format_print=='latex':
                                                        print monomial.latex(),',',
                                                elif format_print=='math':
                                                        print 'Format not implemented'
                                                else:
                                                        print monomial,',',
                                        print '\n','Monomials divisor (general element):',
                                        for monomial in divisor:
                                                if format_print=='latex':
                                                        print monomial.latex(),',',
                                                elif format_print=='math':
                                                        print 'Format not implemented'
                                                else:
                                                        print monomial,',',
                                        print '\n'
                                        if self.is_semistable(M_gamma):
                                                print 'The pair is strictly semistable'
                                                print 'The potential closed orbit associated to this pair is:'
                                                closed_orbit_list=self.closed_orbit[gamma]
                                                for closed_orbit in closed_orbit_list:
                                                        closed_orbit_variety, closed_orbit_divisor=closed_orbit
                                                        print 'N^0(',gamma,')'
                                                        print 'Monomials variety (potential closed orbit):',
                                                        for monomial in closed_orbit_variety:
                                                                if format_print=='latex':
                                                                        print monomial.latex(),',',
                                                                elif format_print=='math':
                                                                        print 'Format not implemented'
                                                                else:
                                                                        print monomial,',',
                                                        print '\n','Monomials divisor (potential closed orbit):',
                                                        for monomial in closed_orbit_divisor:
                                                                if format_print=='latex':
                                                                        print monomial.latex(),',',
                                                                elif format_print=='math':
                                                                        print 'Format not implemented'
                                                                else:
                                                                        print monomial,',',
                                                        print '\n\n'
                                                                
                                        else:
                                                print 'The pair is not semistable'
                                        print '----------------------------------------------\n'
				i=i+1
				print '========================================================\n'

		else:
			for gamma, M_gamma in self.non_stable.iteritems():
				print '(',i,') ',
				print 'N^+(',gamma,')={',
				for monomial in M_gamma:
					print monomial,','
				print '}.'
				print '\n'
				








