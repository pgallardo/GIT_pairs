##################################################################
#Package name: VGIT.Monomials
#Creation date: 15/05/20
#Description: Contains the subpackage Monomials with the class
#       Monomials to deal with one monomial and the class MonFrozenSet to deal with
#       set of monomials.
#Warning: No other subpackage should work directly with the attributes
#       of monomials.
##################################################################
from sympy.matrices import *
from sympy import symbols
from sympy import solve_linear_system
from sympy import Rational
from sympy import *
import itertools
########################
#GLOBAL VARIABLES
########################

#empty

########################
#AUXILIARY FUNCTIONS
########################
def Hilbert_Mumford(pair, gamma, t):
    var,bd=pair
    minimum_var=None
    minimum_bd=None
    for x in var:
        value=gamma.Prod(x)
        if (minimum_var==None):
            minimum_var=value
        else:
            minimum_var=min(minimum_var,value)
    for x in bd:
        value=gamma.Prod(x)
        if (minimum_bd==None):
            minimum_bd=value
        else:
            minimum_bd=min(minimum_bd,value)
    return minimum_var+t*minimum_bd

def Jlcm(*args):
    from sympy.core.numbers import igcd
    """Computes integer least common multiple.

    Examples
    ========

    >>> from sympy.core.numbers import ilcm
    >>> ilcm(5, 10)
    10
    >>> ilcm(7, 3)
    21
    >>> ilcm(5, 10, 15)
    30

    """
    if 0 in args:
        return 0
    a=args[0][0]
    for b in args[0][1:]:
        a = a*b // igcd(a, b)
    return a

def construct_boundary(dim):
    boundary_list=[]
    for i in range(1, dim):
        b=[0 for k in range(dim)]
        b[i]=1
        b[i-1]=-1
        mon=Monomial(b)
        boundary_list.insert(i,mon)
    return boundary_list

def associated1ps(monomials,boundary=None, dim=4): #generates the normalized 1-PS from a list of monomials.
        monomials_copy=list(monomials)
        if(len(monomials_copy)>0):
            m0=monomials_copy[0]
            monomial_list=[m0-m for m in monomials_copy]
        else:
            monomial_list=[]
        #observe that the previous line also adds the zero vector, but it does not change the code.
        
        #we add the monomial in the centre, which is equivalent to saying that lambda is in SL(n,C).
        monomial_list.append(Monomial([1 for x in range(dim)]))

        #We add the boundary components:
        if boundary!=None:
            for x in boundary:
                monomial_list.append(x)
        n_mon=len(monomial_list);
        #building the general matrix
        #first we put the first component with negative value and reverse the monomials
        #to assure that this one becomes the solution vector b, i.e. Ax=-y[0]. This is equivalent to
        #making sure that solution_1=1
        mon_reversed=[]
        for x in monomial_list:
            y=list(x)
            y[0]=-y[0]
            y.reverse()
            mon_reversed=mon_reversed+y
        extended=Matrix(n_mon, dim, mon_reversed)
        variables=list(symbols('x0:%d'%n_mon))

        #finding solutions
        system_solution=solve_linear_system(extended, *variables,  rational=True)
        if system_solution==None:
            return None
        #finding least common denominator
        val=system_solution.values()
        for x in val:
            if not isinstance(x,Rational):
                return None
        denominators=[x.q for x in val if x.p!=0]
        lcm=Jlcm(denominators)
        if not lcm:
            return None
        #clearing solutions
        solution=[int(x*lcm) for x in val]
        solution.append(lcm)
        solution.reverse()
        return solution


########################
#Function: AllMon
#Attributes:
#       dim: dimension of the monomials
#       deg: degree of the monomials
#       gamma (default None): 1-parameter subgroup.
#       allMonomials: list of all Monomials.
#Description: returns a list of all monomials of a given dimension and degree.
#       If gamma!=None, then the list contains only those monomials which are
#       not stable (unstable or semi-stable) with respect to gamma.
#       If gamma=None, then it returns all monomials of the given degree and dimension.
#       It uses the recursive function AllMon_rec (necessary for arbitrary n,d)
########################
def AllMon_rec(dim,deg, gamma,allMonomials, monomial_tup=()):
        if dim>1:
                for i in range(deg+1):
                        AllMon_rec(dim-1, deg-i, gamma, allMonomials, monomial_tup+(i,))
        else:
                add_monomial=Monomial(list(monomial_tup) + [deg])
###*** TO DO: The gamma part has to go in AllMon in order to put the t and both monomials in the product.
                allMonomials.append(add_monomial)
                    
def AllMon(dim,deg, gamma=None, pairs=0, div_degree=1, t=0):
        allMonomials=list()
        allDivisors=list()
        allPairs=list()
        AllMon_rec(dim,deg,gamma,allMonomials)
        if pairs:
            AllMon_rec(dim,div_degree,gamma,allDivisors)
            for pair in itertools.product(allMonomials, allDivisors):
                if not gamma:
                    allPairs.append(pair)
                elif gamma.Prod(pair,pairs, t)>=0:
                    allPairs.append(pair)
            return allPairs
        else:
            if not gamma:
                return allMonomials
            else:
                for monomial in allMonomials:
                    if gamma.Prod(monomial)>=0:
                        allPairs.append(monomial)
                return allPairs

########################
#CLASES
########################

########################
#Class: Monomial
#Attributes:
#       None. Monomial is mainly just a wrap of the class tuple.
#Description: Class of monomials of arbitrary finite dimension.
#Methods:
#   __init__(self,data): creates a Monomial from data
#   __ge__(A,B): redefines the simbol ">=" to be  the partial order induced by 1-PS
#   __le__(A,B): redefines the simbol "<=" to be  the partial order induced by 1-PS
#   deg(self): returns the degree of self
#   Com(self,m2): returns 1 if self>=m_2, -1 if self<=m_2,
#               0 if they are uncomparable
########################


class Monomial(tuple):
        def __init__(self, data):
                self=tuple(data)
        def __hash__(self):
                if len(self)==0:
                        return 0
                else:
                        return self[0]
                
        def __ge__(A,B):
                a_0=b_0=0
                for a_i,b_i in zip(A,B):
                        a_0=a_0+a_i
                        b_0=b_0+b_i
                        if a_0<b_0:
                                break
                else:
                        return True
                return False
        
        def __le__(A,B):
                a_0=b_0=0
                for a_i,b_i in zip(A,B):
                        a_0=a_0+a_i
                        b_0=b_0+b_i
                        if a_0>b_0:
                                break
                else:
                        return True
                return False
        def __sub__(A,B):
                if A.dim()!=B.dim():
                        return None
                output=[a_i - b_i for a_i, b_i in zip(A, B)]
                return output
            
        def deg(self): 
                #reduce adds the first two elements of self, the result to
                #the next two, etc.
                return reduce( lambda x, y: x+y,self)
        def loc_deg(self): 
                return self.deg()-self[len(self)-1]
            
        def Com(self,m2): 
                if (not self>=m2) and (not m2>=self): # true only if not comparable
                        return 0        
                elif (self>=m2):        
                        return 1
                elif (m2>=self):
                        return -1
        def dim(self):
                return len(self)
        def larger_equal(self, MonomialsSet):
                monomialsList=[]
                for monomial in MonomialsSet:
                    if monomial>=self:
                        monomialsList.append(monomial)
                return tuple(monomialsList)
        def latex(self):
            f=1;
            n=len(self)
            variables=list(symbols('x_0:%d'%n))
            for i in range(0,n):
                f=f*variables[i]**self[i]
            return latex(f)
            

########################
#Class: MonFrozenSet
#Attributes:
#       None. MonFrozenSet is just a subclass of the class set.
#Description: Class of  sets of monomials. Used for different purposes.
########################
class MonFrozenSet(frozenset):
    def __new__(cls, dim, deg, data=[], gamma=None, pairs=0, div_degree=1,t=0):
                #if there is data, it creates the set from data.
                if data:
                    return super(MonFrozenSet,cls).__new__(cls,data)
                #if there is no data, then it uses parameters dim and deg
                #the first thing to check is that the set is not meant to be empty
                if dim<=0 or deg<0:
                    return super(MonFrozenSet,cls).__new__(cls,data)
                #If there are indeed monomials, then we create a list of them, and then a set
                #as required by the set construction. We have two options for this. Either
                #a 1PS was passed as a parameter and then we use tkhis to construct the set, or
                #otherwise, we just create a set with all Monomials for a given dimension
                data=AllMon(dim,deg, gamma, pairs, div_degree,t);
                return super(MonFrozenSet,cls).__new__(cls,data)
                
    def __hash__(self):
                return len(self)
    def multiplicity_point(self):
        return min([monomial.loc_deg() for monomial in self])
    def clear_denominators(self, j,p,q,dim):
        output=[]
        for monomial in self:
            adapted=[degree*q*dim for degree in monomial]
            adapted[j]=adapted[j]+p*dim
            output.append(adapted)
        return output
    def minimum(self, dim):
        minimum=0
        for monomial in self:
            for i in range(dim):
                if monomial[i]==1 and i>minimum: #Note that the order is correct since lambda_i>lambda_{i+1}
                    minimum=i
        return minimum


    def degree(self):
        if len(self)==0:
            return -1
        for x in self:
            return x.deg()

