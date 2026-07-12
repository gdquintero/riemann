"""
Distancia optima d_n(alpha)^2 (CORREGIDA: coeficientes reales / norma H^2 genuina).

La distancia en H^2(C_alpha) es  ||F||^2 = (1/2pi) int_{-inf}^inf |F(a+it)|^2 dt.
Para coeficientes REALES c_k, F = sum c_k G_k - E tiene simetria de reflexion y
(1/2pi)int_{-inf}^inf = (1/pi)int_0^inf.  El problema correcto es

    d_n^2 = min_{c real} || sum_{k=2}^n c_k G_k - E ||^2
          = ||E||^2 - btil^T Atil^{-1} btil,
    Atil = Re<G_j,G_k> = Re(M^H M),   btil = Re<E,G_k> = Re(M^H y).

(Minimizar con c COMPLEJO sobre el semieje [0,inf) subestima d_n^2: ese fue el
error de la version anterior.)  s_n usa los coeficientes de Mobius (ya reales).
"""
import numpy as np
from scipy.special import zeta
from prime_tools import mobius_function

C_REF = 0.0461  # sum_rho 1/|rho|^2 = 2 + gamma - log(4 pi)

def build_grid(T_fine, h_fine, T_max, n_geo):
    t1 = np.arange(0.0, T_fine, h_fine)
    t2 = np.geomspace(T_fine, T_max, n_geo)
    t = np.unique(np.concatenate([t1, t2]))
    w = np.zeros_like(t)
    dt = np.diff(t)
    w[1:-1] = 0.5 * (dt[1:] + dt[:-1])
    w[0] = 0.5 * dt[0]
    w[-1] = 0.5 * dt[-1]
    return t, w

def assemble(alpha, t, w, nmax):
    s = alpha + 1j * t
    Z = zeta(s) / s
    ks = np.arange(2, nmax + 1)
    kpow = np.exp(-np.outer(np.log(ks), s)).T
    kinv = (1.0 / ks)[None, :]
    G = (kpow - kinv) * Z[:, None]
    E = 1.0 / s
    sw = np.sqrt(w / np.pi)
    M = G * sw[:, None]
    y = E * sw
    return M, y, ks

def dn2_real(A, b, normE2, n):
    """d_n^2 con Gram real, submatriz k=2..n (indices 0..n-2)."""
    j = n - 1
    Aj = A[:j, :j]
    bj = b[:j]
    return normE2 - float(bj @ np.linalg.solve(Aj, bj))

def compute(alpha, ns, grid_params):
    t, w = build_grid(*grid_params)
    normE2 = np.sum((w / np.pi) * np.abs(1.0 / (alpha + 1j * t)) ** 2)
    M, y, ks = assemble(alpha, t, w, max(ns))
    A = np.real(M.conj().T @ M)          # Gram real (correcta)
    b = np.real(M.conj().T @ y)
    mu = np.array([mobius_function(k) for k in ks], dtype=float)
    dn2, sn = [], []
    for n in ns:
        j = n - 1
        dn2.append(dn2_real(A, b, normE2, n))
        # s_n = || sum mu(k) G_k - E ||^2 = mu^T A mu - 2 mu^T b + ||E||^2
        m = mu[:j]
        sn.append(float(m @ A[:j, :j] @ m - 2.0 * m @ b[:j] + normE2))
    return np.array(dn2), np.array(sn), len(t), normE2

def fit_C(ns, dn2):
    x = 1.0 / np.log(ns)
    Amat = np.vstack([x, x**2]).T
    coef, *_ = np.linalg.lstsq(Amat, dn2, rcond=None)
    return coef

def main():
    alpha = 0.5
    ns = np.array([10, 20, 40, 80, 160, 320, 480, 640, 960, 1280, 1920, 2560, 3200, 4096, 5120])
    grid = (80.0, 0.010, 3.0e5, 20000)

    dn2, sn, gpts, normE2 = compute(alpha, ns, grid)
    C, D = fit_C(ns, dn2)

    print(f"grid points = {gpts}   ||E||^2 = {normE2:.6f} (exacto 1.0)")
    print(f"{'n':>6} {'d_n^2':>13} {'s_n':>13} {'d_n^2*log n':>13} {'s_n*log n':>11}")
    for i, n in enumerate(ns):
        print(f"{n:>6} {dn2[i]:>13.5e} {sn[i]:>13.5e} {dn2[i]*np.log(n):>13.5f} {sn[i]*np.log(n):>11.4f}")
    print()
    print(f"Ajuste d_n^2 ~ C/log n + D/(log n)^2 :  C = {C:.5f}, D = {D:.5f}")
    print(f"Referencia sum_rho 1/|rho|^2 ~ {C_REF}")

    np.savez("dn_results.npz", alpha=alpha, ns=ns, dn2=dn2, sn=sn,
             C=C, D=D, C_ref=C_REF, grid_pts=gpts)
    with open("dn_results.csv", "w") as f:
        f.write("n,dn2,sn,dn2_logn\n")
        for i, n in enumerate(ns):
            f.write(f"{n},{dn2[i]:.8e},{sn[i]:.8e},{dn2[i]*np.log(n):.8e}\n")
    print("Guardado: dn_results.npz, dn_results.csv")

if __name__ == "__main__":
    main()
