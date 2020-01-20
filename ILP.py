
# coding: utf-8

# In[ ]:


from gurobipy import *

W=4     # the width of NoC
H=4     # the hieght of Noc
S=W*H   # the number of core
C=1024  # capacity of a core
N=38    # the number of neurons
M=152   # the number of edges between neurons
x=[]    # the horizontal coordinate of core_i
y=[]    # the vertical coordinate of core_j
for i in range(W):
    for j in range(H):
        x.append(i)
        y.append(j)
l=[0]*N
d=[[0]*N for i in range(N)]
# read in
file=open("computation_graph.txt")
#file=open("graph.txt")
name={}
for i in range(N):
    line=file.readline().strip()
    neuron,str=line.split()
    l[i]=int(str)
    name[neuron]=i
for k in range(M):
    line=file.readline().strip()
    x1,x2,str=line.split()
    i=name[x1]
    j=name[x2]
    d[i][j]=int(str)
    print(i,j,str)
    
try:
    #create a new model
    m=Model("SNN Mapping")
    
    # create variables
    a=m.addVars(N,S,vtype=GRB.BINARY,name="a")  # if neuron_i is mapped into core_j
    f=m.addVars(N,vtype=GRB.INTEGER,name="f")   # the time when neuroni is finished
    b=m.addVars(N,N,vtype=GRB.BINARY,name="b")  # auxiliary variables for the serial core constraint
    t=m.addVars(N,N,vtype=GRB.INTEGER,name="t") # data transfer latency from neuron_i to neuron_j
    ans=m.addVar(vtype=GRB.INTEGER,name="ans")      # maximum f_i

    # set objective
    m.setObjective(ans,GRB.MINIMIZE)

    # add constraint
    m.addConstr(ans==max_(f[i] for i in range(N)),name="ans=maximum f_i")
    #m.addConstr(ans==t.sum('*','*'))
    m.addConstrs((f[i]>=l[i] for i in range(N)),name="finish time is none-negative")
    m.addConstrs((t[i,j]>=0 for i in range(N) for j in range(N)),
                 name="transfer latency is none-negative")
    
    m.addConstrs((a.sum('*',j)<=C for j in range(S)),name="core capacity constraint")
    m.addConstrs((a.sum(i,'*')==1 for i in range(N)),name="no copy constraint")
    m.addConstrs((f[i]+t[i,j]+l[j]<=f[j] for i in range(N) for j in range(N) if d[i][j]>0),
                 name = "neuron order constraint")
    m.addConstrs((-100000*b[i,j]-1000000*(1-a[i,k])-1000000*(1-a[j,k])<=f[i]-l[i]-f[j] 
                  for i in range(N) for j in range(i) for k in range(S)),
                 name = "serial core constraint1")
    m.addConstrs((100000*(b[i,j]-1)-1000000*(1-a[i,k])-1000000*(1-a[j,k])<=f[j]-l[j]-f[i] 
                  for i in range(N) for j in range(i) for k in range(S)),
                 name = "serial core constraint2")
    for i in range(N):
        for j in range(N):
            if d[i][j]<=0:
                m.addConstr(t[i,j]==0)
            expr1=LinExpr(0)
            expr2=LinExpr(0)
            expr3=LinExpr(0)
            expr4=LinExpr(0)
            for k in range(S):
                expr1+=x[k]*a[i,k]-x[k]*a[j,k]+y[k]*a[i,k]-y[k]*a[j,k]
                expr2+=x[k]*a[j,k]-x[k]*a[i,k]+y[k]*a[i,k]-y[k]*a[j,k]
                expr3+=x[k]*a[i,k]-x[k]*a[j,k]+y[k]*a[j,k]-y[k]*a[i,k]
                expr4+=x[k]*a[j,k]-x[k]*a[i,k]+y[k]*a[j,k]-y[k]*a[i,k]
            m.addConstr(d[i][j]*expr1<=t[i,j])
            m.addConstr(d[i][j]*expr2<=t[i,j])
            m.addConstr(d[i][j]*expr3<=t[i,j])
            m.addConstr(d[i][j]*expr4<=t[i,j])
    m.optimize()
    
    for v in m.getVars():
        if v.x>0:
            print(v.varName, v.x)
    print('answer=',m.objVal)
    
except GurobiError:
    print(GurobiError.errno)
    print(GurobiError.message)

