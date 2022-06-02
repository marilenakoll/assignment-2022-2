import argparse

def Compare(A,B,m,d):
    if A==B:
        return m
    else:
        return -d
    

def reverse(A): 
    return A[::-1]

def swap(L,K):
    Ltemp=L.copy() 
    L=K.copy()
    K=Ltemp.copy()
    return L,K
        


def create_F(A,B,g,m,d):
    #create F
    g=abs(g)
    m=abs(m)
    d=abs(d)
    F=[[0 for _ in range(len(B)+1)] for j in range(len(A)+1)]
    for i in range(len(B)+1):
        F[0][i]=-g*i
    for i in range(len(A)+1):
        F[i][0]=-g*i
    for i in range(1,len(A)+1):
        for j in range(1,len(B)+1):
            comp=Compare(A[i-1],B[j-1],m,d) # comp>0 if match, comp<0 if mismatch
            F[i][j]=max(F[i-1][j]-g , F[i][j-1]-g,F[i-1][j-1]+comp)
    return F
    

    

def NeedlemanWunsch(A,B,F,W,Z,g,m,d):  
    global WWnw,ZZnw
    #enumerate alignments
    i=len(A)
    j=len(B)
    g=abs(g)
    m=abs(m)
    d=abs(d)
    if i==0 and j==0:
        WWnw.append(W)
        ZZnw.append(Z)
        return WWnw,ZZnw
    if i>0 and j>0:
        comp=Compare(A[i-1],B[j-1],m,d)
        if F[i][j]==F[i-1][j-1]+comp:
            NeedlemanWunsch(A[0:i-1],B[0:j-1],F,A[i-1]+W,B[j-1]+Z,g,m,d)
    if i>0 and F[i][j]==F[i-1][j]-g:
        NeedlemanWunsch(A[0:i-1],B,F,A[i-1]+W,"-"+Z,g,m,d)
    if j>0 and F[i][j]==F[i][j-1]-g:
        NeedlemanWunsch(A,B[0:j-1],F,"-"+W,B[j-1]+Z,g,m,d)
    



def ComputeAlignmentScore(A,B,g,m,d):
    g=abs(g)
    m=abs(m)
    d=abs(d)
    L=["" for _ in range(len(B)+1)] 
    for j in range(len(L)): 
        L[j]=-j*g
    K=["" for _ in range(len(B)+1)] 
    for i in range(1,len(A)+1):
        L,K=swap(L,K)
        L[0]=-i*g
        for j in range(1,len(B)+1):
            md=Compare(A[i-1],B[j-1],m,d) 
            L[j]=max(L[j-1]-g,K[j]-g,K[j-1]+md)
    return L


def UpdateAlignments(WW,ZZ,WWl,WWr,ZZl,ZZr): 
    WWnew=[]
    ZZnew=[]
    for l in WWl:
        for r in WWr:
            WWnew.append(l+r)
    for l in ZZl:
        for r in ZZr:
            ZZnew.append(l+r)
    for i in range(len(WWnew)):
        WWnew[i]="".join(WWnew[i])
        WW.append(WWnew[i])
    for i in range(len(ZZnew)):
        ZZnew[i]="".join(ZZnew[i])
        ZZ.append(ZZnew[i])
    return WW,ZZ

def element_wise_add(Sl,Sr_rev):
    S=[]
    for i in range(len(Sl)):
        S.append(Sl[i]+Sr_rev[i])
    return S

def Hirschberg(A,B,g,m,d,t):
    global WWnw,ZZnw
    
    if len(A)==0:
        WW=["-"*len(B)] 
        ZZ=[B]
        
    elif len(B)==0:
        WW=[A]
        ZZ=["-"*len(A)]
    elif len(A)==1 or len(B)==1:
        F=create_F(A,B,g,m,d)
        NeedlemanWunsch(A,B,F,"","",g,m,d) 
        WW=WWnw 
        ZZ=ZZnw
        WWnw=[]
        ZZnw=[]
    else:
        i=int(len(A)/2) 
        Sl=ComputeAlignmentScore(A[0:i],B,g,m,d)
        Arev=reverse(A[i:len(A)])
        Brev=reverse(B)
        Sr=ComputeAlignmentScore(Arev,Brev,g,m,d) 

        Sr_rev=reverse(Sr)

        S=element_wise_add(Sl,Sr_rev)
        J=max(S) 
        J_list=[]
        for k,score in enumerate(S):
            if score==J:
                J_list.append(k)
        WW=[]
        ZZ=[]
        for j in J_list:
            if args.t:
                print(i,", ",j,sep="")
            WWl,ZZl=Hirschberg(A[0:i],B[0:j],g,m,d,t)
            WWr,ZZr=Hirschberg(A[i:len(A)],B[j:len(B)],g,m,d,t)
            WW,ZZ=UpdateAlignments(WW,ZZ,WWl,WWr,ZZl,ZZr)
    return WW,ZZ




parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store_true')
parser.add_argument('g', type=int)
parser.add_argument('m', type=int)
parser.add_argument('d', type=int)
parser.add_argument('A', type=str)
parser.add_argument('B', type=str)
args = parser.parse_args()


WWnw=[]
ZZnw=[]
WW,ZZ=Hirschberg(args.A,args.B,args.g,args.m,args.d,args.t)
results=[]
WWset=set(WW)
ZZset=set(ZZ)
if len(WWset)!=len(ZZset):
    for i in WW:
        for j in ZZ:
            if [i,j] not in results:
                results.append([i,j])
    for r in results:
        print(r[0])
        print(r[1])
        if r!=results[-1]:
                print("")
else:
    
    for i in range(len(WW)):
        if [WW[i],ZZ[i]] not in results:
            results.append([WW[i],ZZ[i]])
    for r in results:
            print(r[0])
            print(r[1])
            if r!=results[-1]:
                print("")
