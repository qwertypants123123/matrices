import random
from numpy import polynomial as P
import copy, math

class poly:
    def __init__(self,l,name="x"):
        self.l=l#like in numpy, coefficients are stored backwards(ex. 1+2x+3x^2+...)
        self.name=name
        self.le=len(l)
        #for i in range(len(l)): #for if i ever try to use fractions
            #if type(l[i])==int:
                #l[i]=frac(l[i])
    
    def __str__(self):
        if len(self.l) == 0:
            return "0"
        ret=""
        l=list(self.l)
        for i in range(len(l) - 1, -1, -1):
            if l[i] == 0 and i != 0:#omit zero coefficients
                continue
            if l[i] > 0:#join with a plus symbol, if pos. if neg, already included in string version of coeff
                ret += '+'
            ret += (f'{l[i]}{self.name}^{i}')
        if ret[0] == '+': #remove + from first, unnecessary
            ret = ret[1:]
        ret = ret[:len(ret)-3]#from 1x^0, remove x^0
        return ret
        
    def __mul__(a,b):#b has to be the integer multiplier
        if type(a)==int:
            l=b.l
            for i in range(b.le):
                l[i]*=a
            return poly(l,a.name)
        if type(b)==int:
            l=a.l
            for i in range(a.le):
                l[i]*=b
            return poly(l,a.name)
        if type(a)==poly and type(b)==poly:
            ret=[0] * (a.le+b.le-1)#degree of product is equal to sum of degrees of a, b. 
            for i in range(a.le):
                for i2 in range(b.le):
                    ret[i+i2]+=a.l[i]*b.l[i2]#starts from x^0 on both sides, so product is x^(a+b)
            return poly(ret,a.name)
    
    def __add__(a,b):
        if type(a)==poly and type(b)==poly:
            if len(a.l)>len(b.l):
                l=list(a.l)
                li=list(b.l)
            else:
                return b + a
            for i in range(len(li)):
                l[i]+=li[i]
            return poly(l,a.name)
        if type(b) == float or type(b) == int:
            a.l[0] += b
            return poly(a.l)
        else:
            return b + a
    
    def __neg__(self):
        l = list(self.l)
        for i in range(len(l)):
            l[i] = -1 * l[i]
        return poly(l)
    
    def __sub__(a, b):
        return a + (-b)

class matrix:
    def __init__(self,l=None, round = True): 
        if l!=None:
            self.special=set([])
            self.le=len(l[0])
            for i in l:
                assert len(i)==self.le,'must be a rectangular matrix'
            self.he=len(l)
            self.l=l
            if self.le==self.he:
                self.special.add("square")
            if self.le==1:
                self.special.add("vector")
            if round:
                for i in range(self.he):
                    for i2 in range(self.le):
                        try:
                            if math.isclose(abs(l[i][i2]), 0, rel_tol=1e-10):#removes issues with some float operations
                                l[i][i2] = 0
                        except:
                            1==1
#                    if l[i][i2]==int:
#                        l[i][i2]=(l[i][i2])

    def __str__(self):
        return f'{self.l}'
    
    def leading_zeroes(l):
        count = 0
        for i in l:
            if math.isclose(abs(i), 0, abs_tol = 1e-10):#ignores floating point weird stuff
                count += 1
            else:
                break
        return count
    
    def __sub__(self, other):
        try:
            assert self.le == other.le and self.he == other.he, 'invalid matrix sizes'
        except:
            raise TypeError#not matrices
        ret = []
        for i in range(self.he):#per row
            temp = []
            for j in range(self.le):#per col
                temp.append(self.l[i][j] - other.l[i][j])
            ret.append(temp)
        return matrix(ret)
    
    def __mul__(a,b):
        if type(b)==int:
            for i in range(b.le):
                for i2 in range(b.he):
                    a.l[i][i2]*=b
        if type(a) == int:
            return b * a
        if type(b)==matrix:
            assert a.le==b.he,'invalid matrix size'
            ret=[]
            for i in range(a.he):#per row of output
                temp=[]
                for i2 in range(b.le):#per column of output
                    sum=0
                    for i3 in range(a.le):#per index in multiplication
                        sum+=a.l[i][i3]*b.l[i3][i2]
                    temp.append(sum)
                ret.append(temp)
            return ret
        
    def inverse(self):
        try:
            assert "square" in self.special,'can only take inverses of a square matrix'
            le=self.le
            ide=matrix(matrix.identity(le))
            he=self.he
            l=self.l
            de=l        
            for i in range(1,le):#each row
                for i2 in range(i):#each entry
                    m=l[i][i2]/l[i2][i2]
                    self.subtract(i,i2,m)
                    ide.subtract(i,i2,m)
            for i in range(le):
                p=l[i][i]
                for i2 in range(le):
                    l[i][i2]/=p
                    try:
                        ide.l[i][i2]/=p
                    except:
                        print(ide.l)
            
            for i in range(le-2,-1,-1):
                for i2 in range(le-1,i,-1):
                    e=l[i][i2]
                    ide.subtract(i,i2,e)
            for i in range(ide.he):
                    for i2 in range(ide.le):
                        try:
                            ide.l[i][i2]=ide.l[i][i2].print()
                        except:
                            None
            return matrix(ide.l)
        except ZeroDivisionError:
            return("Cannot take inverse of singular matrix")
    
    def identity(x, multiple = 1):
        l=[]
        for le in range(x):
            l2=[]
            for w in range(x):
                if le==w:
                    l2.append(multiple)
                else:
                    l2.append(0)
            l.append(l2)
        return matrix(l)
    
    def rand_matrix(h, l):
        ret = []
        for i in range(h):
            temp = []
            for j in range(l):
                temp.append(random.random())
            ret.append(temp)
        return matrix(ret)

    def subtract(self, f, l, mult):#subtracts a multiple of one row from another
        le=self.le
        mainl=self.l
        for i in range(le):
            mainl[f][i]-=(mainl[l][i]*mult)
    
    def transpose(self):
        ret=[]
        for i in range(self.le):
            temp=[]
            for i2 in range(self.he):
                temp.append(self.l[i2][i])
            ret.append(temp)
        return ret

    def det(self):
        le=self.le
        he=self.he
        l=list(self.l)
        vals=[]
        temp=[]
        cofactors=[]
        sum=poly([0])
        temp2=[]
        if he==2:
            return (l[0][0]*l[1][1])-(l[0][1]*l[1][0])
        else:
            for i in range(le):
                vals.append(i)
            for i in vals:
                for i2 in range(1,le):# per row, so one less than height of matrix
                    for i3 in range(le):# per column, same but one is elim thr next line
                        if i3!=i: # if cofactor thing is not on same column as entry, append
                            temp.append(l[i2][i3])
                    temp2.append(temp)
                    temp=[]
                exec(f'cofactors.append(matrix(temp2))')
                temp2=[]
            for i in range(le):
                sum+=(l[0][i]*cofactors[i].det()*((-1)**(i)))#multiplying by 1, -1, and by cofactor thing
            return sum
    def eigvals(self):
        temp=copy.deepcopy(self.l)#required because list contains mutable objects(lists), not just copying references
        assert "square" in self.special
        for i in range(self.he):
            for i2 in range(self.le):
                if i == i2:
                    temp[i][i2] = poly([temp[i][i2], -1])
                else:
                    temp[i][i2] = poly([temp[i][i2]])
        temp=matrix(temp)
        temp = temp.det()
        p=P.Polynomial(temp.l)
        return p.roots()
    
    def ref(self):
        l = self.l
        cont = True
        height = self.he
        col = 0
        row = 0 #pivot numbers
        l = sorted(l, key = matrix.leading_zeroes)
        self.l = l
        while cont:
            for i in range(row + 1, height):
                mult = l[i][col] / l[row][col]
                self.subtract(i, row, mult)
                for j in range(self.le):
                    if math.isclose(abs(l[i][j]), 0, abs_tol = 1e-10):
                        l[i][j] = 0
            rem = l[row + 1:]#remaining rows
            rem = sorted(rem, key = matrix.leading_zeroes, reverse = True)
            row += 1
            zeroes = matrix.leading_zeroes(l[row])
            if zeroes != self.le:
                col = zeroes
            if row == self.he - 1 or zeroes == self.le:
                cont = False
        return self
    
    def rref(self):
        self.ref()
        l = self.l
        height = self.he
        row, col = 0, 0
        for i in range(height):
            if l[row][col] != 0:
                for j in range(0, row):
                    mult = l[j][col] / l[row][col]
                    self.subtract(j, row, mult)
                mult = 1 / l[row][col]
                for j in range(col, self.le):
                    l[row][j] *= mult
                l[row][col] = 1
                row += 1
            col += 1
        return self
    
    def nullspace(self):
        self.rref()
        l = self.l
        m = self.le
        n = self.he
        zeroes = [0] * n
        pivots = []
        ret = []
        row = 0
        while row < n:#appends column numbers for pivots
            if row != 0:
                for i in range(row, n):
                    if l[row][i] == 1 and l[row - 1][i] == 0:
                        pivots.append(i)
                        break
            else:
                for i in range(row, n):
                    if l[row][i] == 1:
                        pivots.append(i)
                        break
            row += 1
        for i in range(n):#going over columns
            if i not in pivots:#for free variables:
                temp = [0] * m
                for j in range(len(pivots)):#rows
                    temp[pivots[j]] -= l[j][i]#given a positive one times a free variable, neg the others
                if temp != zeroes:
                    ret.append(temp)
                    ret[len(ret) - 1][i] += 1
                else:
                    ret.append(list(zeroes))
                    ret[len(ret) - 1][i] += 1
        return ret
        
    def eigvectors(self):#doesn't work with complex eig's
        eigvals = self.eigvals()
        l = []
        width = self.le
        ret = []
        for i in eigvals:
            l.append(complex(i))
        for i in l:
            I = matrix.identity(width, i)
            temp = self - I
            nullspace = temp.nullspace()
            for j in nullspace:
                ret.append(j)
        return ret#success!!
    
    #def __pow__(self, other):
        


a = matrix([
    [0,   2.5, -1.3],
    [-2.5, 0,   4.7],
    [1.3, -4.7, 0 ]
])
print(a.eigvals())