"""
Barrido en alpha del sensor de Theta (Corolario 3.2).

Para cada alpha computamos la distancia optima d_n(alpha)^2 (Gram real, estable) y la
comparamos con la COTA INFERIOR del Teorema 3.2 (piso por ceros):
    d_n(alpha)^2 >= v* K^{-1} v,   K_ij = 1/(rho_i + conj(rho_j) - 2 alpha),  v_i = 1/rho_i,
usando los ceros rho = 1/2 + i gamma con Re(rho) > alpha (i.e. solo si alpha < 1/2).

Idea: para alpha < 1/2 hay ceros con Re > alpha  ->  d_n^2 NO ->0 (piso positivo).
Para alpha >= 1/2 no hay obstruccion  ->  d_n^2 decae. La transicion esta en alpha=1/2.
Guarda dn_sweep.npz.
"""
import numpy as np
from scipy.special import zeta

# partes imaginarias de los primeros ceros no triviales de zeta
GAMMAS = np.array([
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
    67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
])

def build_grid(T_fine=80.0, h_fine=0.01, T_max=3.0e5, n_geo=15000):
    t1 = np.arange(0.0, T_fine, h_fine)
    t2 = np.geomspace(T_fine, T_max, n_geo)
    t = np.unique(np.concatenate([t1, t2]))
    w = np.zeros_like(t)
    dt = np.diff(t)
    w[1:-1] = 0.5 * (dt[1:] + dt[:-1])
    w[0] = 0.5 * dt[0]
    w[-1] = 0.5 * dt[-1]
    return t, w

def gram(alpha, t, w, nmax):
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
    A = np.real(M.conj().T @ M)
    b = np.real(M.conj().T @ y)
    normE2 = float(np.real(np.vdot(y, y)))
    return A, b, normE2

def dn2_curve(A, b, normE2, ns):
    return np.array([normE2 - float(b[:n-1] @ np.linalg.solve(A[:n-1, :n-1], b[:n-1]))
                     for n in ns])

def floor_bound(alpha, mzeros=len(GAMMAS)):
    """Cota inferior del Teorema 3.2 usando ceros con Re(rho)=1/2 > alpha."""
    if alpha >= 0.5:
        return 0.0
    rho = 0.5 + 1j * GAMMAS[:mzeros]
    K = 1.0 / (rho[:, None] + np.conjugate(rho)[None, :] - 2.0 * alpha)
    v = 1.0 / rho
    val = np.real(np.conjugate(v) @ np.linalg.solve(K, v))
    return float(val)

def main():
    alphas = np.array([0.35, 0.40, 0.44, 0.47, 0.50, 0.53, 0.56, 0.60,
                       0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.98])
    ns = np.array([20, 40, 80, 160, 320, 640, 1280, 2560])
    t, w = build_grid()
    print(f"grid points = {len(t)}")

    DN2 = np.zeros((len(alphas), len(ns)))
    floor = np.array([floor_bound(a) for a in alphas])
    for i, a in enumerate(alphas):
        A, b, nE2 = gram(a, t, w, ns.max())
        DN2[i] = dn2_curve(A, b, nE2, ns)
        # tasa aparente log-log en la cola
        m = ns >= 160
        gam = -np.polyfit(np.log(ns[m]), np.log(DN2[i][m]), 1)[0]
        print(f"alpha={a:.2f}  ||E||^2={nE2:.4f}(exact {1/(2*a):.4f})  "
              f"d_2560^2={DN2[i,-1]:.3e}  piso(Thm3.2)={floor[i]:.3e}  gamma~{gam:.3f}")

    np.savez("dn_sweep.npz", alphas=alphas, ns=ns, DN2=DN2, floor=floor)
    print("Guardado: dn_sweep.npz")

if __name__ == "__main__":
    main()
