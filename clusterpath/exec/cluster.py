from math import sqrt
from cvxmod import matrix,concatvert,problem,optvar,rows,minimize,concathoriz
from cvxmod.atoms import norm1,norm2,norminf
def rtab(t):
    return zip(*[[float(x) for x in l.split()] for l in open(t).readlines()])
def n1(V):
    return sum([abs(x) for x in V])
def sum_squares(V):
    return sum([x**2 for x in V])
def n2(V):
    return sqrt(sum_squares(V))
def ninf(V):
    return max([abs(x) for x in V])
def getQij(qij,P):
    """Construct the matrix used for diff bt alpha i and j.

    We return Qij and then do sum_{i\neq j} norm(Qij*alpha)

    """
    Qij=matrix()
    zero=matrix(0.0,qij.size)
    for k in range(P):
        thisrow=matrix()
        for l in range(P):
            new = qij if k==l else zero
            thisrow = concathoriz(thisrow,new)
        Qij=concatvert(Qij,thisrow)
    return Qij
def solve(cols,W,penalty,par,parvalue):
    P=len(cols)
    N=len(cols[0])
    norm=eval("norm"+penalty)
    normFUN=eval("n"+penalty)
    normQX=0.0
    SUMDIFF=0.0
    X=matrix([list(x) for x in cols],(N*P,1))
    alpha=optvar("alpha",N*P)
    PARAM=float(parvalue)
    for i in range(N):
        for j in range(i+1,N):
            w=W[i][j]
            q=matrix(0.0,(1,N))
            q[i]=1
            q[j]=-1
            Qij=getQij(q,P)
            SUMDIFF+= norm(Qij*alpha)*w
            xidiff=[cols[k][i]-cols[k][j] for k in range(P)]
            normQX += normFUN(xidiff)*w
    from cvxmod.atoms import sum,square
    ## TDH alternate parameterization
    ##penalty_norm = 1.0/float(N*(N-1))
    ##SUMDIFF *= penalty_norm
    ##normQX  *= penalty_norm
    ##error=(0.5/N)*sum(square(X-alpha))
    error=0.5*sum(square(X-alpha))
    problems={
        "lambda":(error+PARAM*SUMDIFF,[]),
        "s":(error,[SUMDIFF*(1/normQX) <= PARAM]),
        }
    tomin,constraint=problems[par]
    p=problem(minimize(tomin),constraint)
    p.solve()
    return alpha

if __name__ == "__main__":
    import sys
    try:
        table=sys.argv[1]
        par,parvalue=sys.argv[2].split("=")
        penalty=sys.argv[3]
        outfile=sys.argv[4]
    except:
        print "Usage: %s table (s|lambda)=value (1|2|inf) out"%sys.argv[0]
        print sys.argv
        sys.exit(1)
    cols=rtab(table)
    W=rtab(table+".W")
    alpha=solve(cols,W,penalty,par,parvalue)
    o=open(outfile,"w")
    o.write('\n'.join([str(x) for x in alpha.value]))
