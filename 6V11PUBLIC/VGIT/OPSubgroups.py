##################################################################
#Package name: VGIT.OPSubgroups
#Creation date: 15/05/20
#Description: Contains the subpackage 1PSubgroups with the class
#       1PS and all the methods to work with them.
#Warning: No other subpackage should work directly with the attributes
#       of 1PS.
##################################################################
import itertools
from sympy import Rational
from Monomials import associated1ps
from Monomials import MonFrozenSet
########################
#GLOBAL VARIABLES
########################
hash_constant=1000
#empty

########################
##AUXILIARY FUNCTIONS
## returns 1 if semistable, 0 if unstable
########################
def semistable(dim, deg, t, M_gamma):
    from copy import deepcopy
    from scipy.spatial import ConvexHull
    import numpy as np
    empty = np.empty([0,dim], dtype=np.dtype(int))
    variety, divisor=M_gamma
    j=divisor.minimum(dim)
    M_gamma_adapted=variety.clear_denominators(j,t.p, t.q, dim)
    debug_centroid=[t.q*deg+t.p for i in range(dim)]
    aid=[0 for i in range(dim)]
    aid[dim-1]=2*deg
    roof=[aid]
    for i in range(dim-1):
        temp=[0 for j in range(dim)]
        temp[i]=2*deg
        roof.append(temp)
    empty=np.append(empty, roof, axis=0)
    points=np.append(empty,M_gamma_adapted, axis=0)
    points2=np.append(points, [[t.q*deg+t.p for i in range(dim)]], axis=0)
    hull = ConvexHull(points) #The option is to allow situations where the convex hull is not full dimensional
    vertices=tuple(tuple(x) for x in points[hull.vertices].tolist())
    hull2= ConvexHull(points2) #The option is to allow situations where the convex hull is not full dimensional
    vertices2=tuple(tuple(x) for x in points2[hull2.vertices].tolist())
    vertices_set=set(vertices)
    vertices_set2=set(vertices2)
    if(vertices_set==vertices_set2):
        return 1
    else:
        return 0



def OPS_destabilized_pair(dim, deg, gamma, t, allMonomials, allDivisors, support_divisor):
    divisor_list=[]
    variety_list=[]
    for element in allDivisors:
        if element>=support_divisor:
            divisor_list.append(element)
    for element in allMonomials:
        if gamma.Prod_pair(element, support_divisor, t)>=0:
            variety_list.append(element)
    return (MonFrozenSet(dim, deg, variety_list), MonFrozenSet(dim, deg, divisor_list))
    


########################
#CLASES
########################

########################
#Class: OPS
#Attributes:
#       None. 1PS is just a wrap of the class tuple.
#Description: Class of diagonal 1-parameter subgroups in SL(n,C) (good enough for
#       GIT problems by Hilbert-Mumford Numerical Criterion, see Mukai Theorem 3.7)
#       Some methods may force that if t is 1PS, then t[i]>=t[i+1] for all i and that
#       sum(t[0],...,t[n-1])=0.
#Methods:
#   __init__(self,data): creates a 1PS from data, assuming
#               sum(data[0],...,data[n-1])==0, else returns None
#   __eq__(A,B): A==B, None if A==0 or B==0
#########################


class OPS(tuple): 
        def __new__(self, data, normalized=True, SL_N=True):
                #If it is not in SL(n,C)  it returns None
                if not data:
                        return None
                if (SL_N):
                        if (reduce(lambda x, y: x+y,data)!=0):
                                self=None
                                return None
                #If the entries are not ordered, it returns None:
                if(normalized):
                        for i in range(len(data)-1):
                                if data[i]<data[i+1]:
                                        self=None
                                        return None
                self = super(OPS, self).__new__(self, data)
                return self
        def __hash__(self):
                if len(self)==0:
                        return 0
                if self[0]==0:
                        return 0
                return self[1]*hash_constant//self[0]

	def __eq__(A,B):
		# multiplying the entries in the object OPS by a constant keeps the 1-PS in the same
		# way because the 1ps is a subgroup of SL(n,C). 
		# Therefore we check equality up to multiplication by a constant,  e.g. A._a/A._b=B._a/B._b iff A.a*B._b=A._b*B._a
                if len(A)!=len(B):
                        return False
#               The following makes a1 to be 0 if and only if a1 is the identity
                a1=reduce(lambda x, y: x or y,A);
                b1=reduce(lambda x, y: x or y,B);
                if a1==0 and b1==0:
                        return True;
                elif a1==0 or b1==0:
                        return False;
                for i in range(len(A)):
                        for j in range(len(B)):
                                if (A[i]*B[j]-A[j]*B[i]!=0):
                                        return False
                return True
		
	def is_Normalized(self):
                if reduce(lambda x, y: x+y,self)!=0:
                        return False
                for i in range(len(self)-1):
                        if self[i]<self[i+1]:
                                return False
                else:	
	    		return True									
		
        def Dual(self): # returns dual of 1-PS L
                return OPS(tuple(-i for i in reversed(self)))



 
	def Prod_pair (self, variety, divisor, t=0):
                return self.Prod((variety, divisor), 1, t)

	def Prod(self,monomial,pairs=0,t=0): # pairing 1PS and monomial
                result=0
                if not pairs:
                        for lambda_i,variety_i in zip(self,monomial):
                                result=result+lambda_i*variety_i
                else:
                        variety, divisor=monomial
                        for lambda_i,variety_i,divisor_i in zip(self,variety,divisor):
                                result=result+lambda_i*(t.q*variety_i + t.p*divisor_i)
                return result
        def Ann(self, pairs, t=0, minimum_hit=0):
                Ann_list=[]
                for pair in pairs:
                        if self.Prod(pair, 1, t)==minimum_hit:
                                Ann_list.append(pair)
                return tuple(Ann_list);
        def dimension(self):
                return len(self);
        
        def annihilator(self, dim, deg, M_gamma_pair, t):
                variety, divisor=M_gamma_pair
                annihilator_variety=set()
                annihilator_divisor=set()
                for monomial in divisor:
                    flag=1
                    for monomial_test in divisor:
                        if self.Prod(monomial_test)<self.Prod(monomial):
                            flag=0
                    if flag==1:
                        annihilator_divisor.add(monomial)
                for  v, b in itertools.product(variety, divisor):
                    if self.Prod_pair(v,b,t)==0:
                        annihilator_variety.add(v)
                return (MonFrozenSet(dim, deg, annihilator_variety), MonFrozenSet(dim, deg, annihilator_divisor))

            

	
	
			



########################
#Class: OPS_set
#Attributes:
#       None. OPS_set is just a wrap of the class set.
#Description: Class of sets of normalized 1-parameter subgroups in SL(n,C)
#Methods:
#   __init__(self, dim, monomialSet): It creates a set with all normalized 1-parameter
#               subgroups associated to dim-1 monomials in monomialSet and the center monomial. We assume all such monomials
#               have dim variables.
#########################


class OPS_set(set):
        def __init__(self,dim, varietySet, boundaryOPS, pairs=0, divisorSet=None):# generates list B of all normalized 1-PS associated to two monomials + center on a configuration self
                if ((not varietySet) or (pairs and divisorSet==None)):
                        return None
                data=[]
                for k in range(0,dim-1):
                    test=1
                    for boundarySubset in itertools.combinations(boundaryOPS,k):                        
                        for monomialSubset in itertools.combinations(varietySet,dim-1-k):
                                A=OPS(associated1ps(list(monomialSubset), list(boundarySubset),dim))
                                if A!=None and A.is_Normalized():
                                        data.append(A)
                super(OPS_set,self).__init__(data)
        def print_set(self):
                print 'The list of 1-parameter subgroups in SL(n) which ',
                print 'are normalized (i.e. with descending weights) is:'
                for x in self:
                        print x

class OPS_dictionary(dict):
        def __init__(self, OPS_list, pairs_list, t):
                if not(OPS_list and pairs_list):
                        return None
                super(OPS_dictionary, self).__init__(data)
                for gamma in OPS_list:
                        for pair in pairs_list:
                                if gamma.Prod(pair, t)==0:
                                        self[gamma].add(pair)
                return self

        
