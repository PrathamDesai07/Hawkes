"""Problem 3 (rest): general-p moments of Z and exponential tail.
Claim 1: E Z^p <= C_p eps^-p, C_0=1, C_p=(1/(p b)) sum_{i=1..p} a_{p,i} C_{p-i}, b=beta_min, B=beta_max,
         a_{p,1}=p mu + C(p,2) B,  a_{p,i}=mu C(p,i)+B C(p,i+1) (i>=2).
Claim 2: E exp(eps Z/e) <= 2 exp(4 mu/(e b)).
Exact stationary moments from the closed linear system for E[S^a]."""
import numpy as np, itertools
from math import comb, e, exp

def exact_moments(mu, beta, alpha, p):
    K = len(beta)
    monos = [a for a in itertools.product(range(p+1), repeat=K) if sum(a) <= p]
    idx = {a: i for i, a in enumerate(monos)}; zero = tuple([0]*K)
    unk = [a for a in monos if a != zero]; col = {a: i for i, a in enumerate(unk)}
    A = np.zeros((len(unk), len(unk))); rhs = np.zeros(len(unk))
    for r, a in enumerate(unk):
        def add(d, coef):
            if d == zero: rhs[r] -= coef
            else: A[r, col[d]] += coef
        add(a, -sum(a[k]*beta[k] for k in range(K)))
        for d in itertools.product(*[range(ak+1) for ak in a]):
            if d == a: continue
            coef = np.prod([comb(a[k], d[k]) for k in range(K)])
            add(d, mu*coef)
            for m in range(K):
                d2 = list(d); d2[m] += 1; add(tuple(d2), alpha[m]*coef)
    sol = np.linalg.solve(A, rhs)
    return lambda a: 1.0 if a == zero else sol[col[a]]

def poly_pow_moment(Em, coefs, const, q, K):
    """E (const + sum coefs_k S_k)^q via polynomial expansion."""
    poly = {tuple([0]*K): 1.0}
    for _ in range(q):
        new = {}
        for a, c in poly.items():
            new[a] = new.get(a, 0) + c*const
            for k in range(K):
                a2 = list(a); a2[k] += 1; a2 = tuple(a2); new[a2] = new.get(a2, 0) + c*coefs[k]
        poly = new
    return sum(c*Em(a) for a, c in poly.items())

def Cp(mu, b, B, pmax):
    C = [1.0]
    for p in range(1, pmax+1):
        s = 0.0
        for i in range(1, p+1):
            a = p*mu + comb(p, 2)*B if i == 1 else mu*comb(p, i) + B*(comb(p, i+1) if i+1 <= p else 0)
            s += a*C[p-i]
        C.append(s/(p*b))
    return C

if __name__ == "__main__":
    mu = 1.0; P = 4
    print("m_p = E Z^p exact vs C_p eps^-p  (ratio = m_p eps^p / C_p, must be <= 1)")
    print(f"{'beta':>10} {'eps':>6} " + " ".join(f"{'p='+str(p):>9}" for p in range(1, P+1)))
    for beta in (np.array([1.]), np.array([1., 2.]), np.array([1., 4.]), np.array([1., 2., 4.])):
        K = len(beta); b, B = beta.min(), beta.max(); C = Cp(mu, b, B, P)
        for eps in (0.3, 0.1, 0.03):
            alpha = (1-eps)*beta/K; c = alpha/beta
            Em = exact_moments(mu, beta, alpha, P)
            r = [poly_pow_moment(Em, c, 0.0, p, K)*eps**p/C[p] for p in range(1, P+1)]
            print(f"{str(beta.tolist()):>10} {eps:6.2f} " + " ".join(f"{x:9.4f}" for x in r))
    print("\nC_p (mu=1,b=1,B=1):", [round(x, 2) for x in Cp(1, 1, 1, 6)], "  vs p!(?)", [1, 1, 2, 6, 24, 120, 720])
