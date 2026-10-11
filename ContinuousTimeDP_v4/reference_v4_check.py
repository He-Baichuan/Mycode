"""SciPy-only algebra checks; NOT a MATLAB execution or HJBkit integration test."""
import numpy as np
from scipy.linalg import solve

# Noncontrolled stopping VI (rho=.1, f=2a, stop=10).
a = np.linspace(0, 1, 61)
rho = 0.1
f = 2 * a
G = np.ones_like(a) * 10
V = np.maximum(f/rho, G)
continuation = f - rho*V
R = np.maximum(G-V, continuation)
assert np.max(np.abs(R)) < 1e-12

# Regime switching via impulse: no flow in regime 1, flow=2 in regime 2.
# one-way intervention to regime 2 with immediate cost 3.
W = np.stack([np.full(61, 17.), np.full(61, 20.)], axis=-1)
impulse = W[:,1]-3
R1 = np.maximum(impulse-W[:,0], -rho*W[:,0])
R2 = 2-rho*W[:,1]
assert np.max(np.abs(R1)) < 1e-12
assert np.max(np.abs(R2)) < 1e-12

# Robust stopper: nature chooses .2 from high state to stop state V1=15.
V2 = (2+0.2*15)/(.1+.2)
assert abs(V2-50/3) < 1e-12
F1 = min([h*(V2-15)-.1*15 for h in [0.1,0.2]])
F2 = min([2+h*(15-V2)-.1*V2 for h in [.1,.2]])
assert F1 < 0 and abs(F2) < 1e-12
print('SciPy independent algebra identities PASS (not MATLAB verification).')
