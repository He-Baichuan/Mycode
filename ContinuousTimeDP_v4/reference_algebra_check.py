"""Independent SciPy algebra checks, NOT an execution of MATLAB sources."""
import numpy as np
from scipy.sparse import diags, eye, kron, block_diag, csr_matrix
from scipy.sparse.linalg import spsolve

rho=0.07
Q=np.array([[-0.10,0.10],[0.04,-0.04]])
f=np.array([1.,2.])
expected=np.linalg.solve(rho*np.eye(2)-Q,f)
N=51
a=np.linspace(0,1,N)
h=1/(N-1)
mu=.2*(.5-a)
up=np.maximum(mu,0)/h; down=np.maximum(-mu,0)/h
up[-1]=0; down[0]=0
A=diags([down[1:],-(up+down),up[:-1]],[-1,0,1],shape=(N,N),format='csr')
assert abs(A.sum(axis=1)).max() < 1e-12
As=block_diag([A,A],format='csr')+kron(csr_matrix(Q),eye(N,format='csr'),format='csr')
V=spsolve(rho*eye(2*N,format='csr')-As,np.r_[np.ones(N),2*np.ones(N)]).reshape(2,N)
assert max(abs(V[0]-expected[0])) < 1e-11
assert max(abs(V[1]-expected[1])) < 1e-11
print('Markov switching benchmark:',expected,'max error',np.max(abs(V-expected[:,None])))

rho=.1; lam=.4; target=.5
a=np.linspace(0,1,55); N=len(a)
P=np.tile(np.eye(N)[np.argmin(abs(a-target))],(N,1))
J=lam*(csr_matrix(P)-eye(N,format='csr'))
v=spsolve(rho*eye(N,format='csr')-J,a)
exact=(a+lam*target/rho)/(rho+lam)
assert np.max(abs(v-exact))<1e-11
print('Poisson reset benchmark max error:',np.max(abs(v-exact)))

# Constant flow, diffusion in both coordinates. Positive-type condition.
N1,N2=12,14; h1,h2=1/(N1-1),1/(N2-1)
c11=.04; c22=.04; c12=.2*.03
assert c11 >= abs(c12)*h1/h2 and c22 >= abs(c12)*h2/h1
print('Correlated diffusion: positive-type dominance verified')
