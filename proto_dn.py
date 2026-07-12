"""
Prototipo Sugerencia 1: distancia optima d_n(alpha)^2 vs coeficientes mu(k).
- Valida la cuadratura con ||E||^2_alpha = 1/(2 alpha).
- Calcula d_n^2 = min_c || Sum c_k G_k - E ||^2  (minimos cuadrados en la frontera).
- Compara con s_n = || Sum mu(k) G_k - E ||^2.
- Reporta d_n^2 * log n para ver si se estabiliza (ley c/log n).
"""
import numpy as np
from scipy.special import zeta
from prime_tools import mobius_function

def build_grid(T_fine=60.0, h_fine=0.02, T_max=2.0e5, n_geo=4000):
    # malla fina cerca de 0 + cola geometrica; pesos trapezoidales
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
    Z = zeta(s) / s                      # zeta(s)/s en la frontera
    ks = np.arange(2, nmax + 1)
    # G_k(s) = (k^{-s} - k^{-1}) * zeta(s)/s
    kpow = np.exp(-np.outer(np.log(ks), s)).T      # (grid, ncol) = k^{-s}
    kinv = (1.0 / ks)[None, :]
    G = (kpow - kinv) * Z[:, None]                 # (grid, ncol)
    E = (1.0 / s)                                  # (grid,)
    sw = np.sqrt(w / np.pi)                        # ||F||^2 = (1/pi) int_0^inf |F|^2
    M = G * sw[:, None]
    y = E * sw
    return M, y, ks

def main():
    alpha = 0.5
    t, w = build_grid()
    normE2 = np.sum((w / np.pi) * np.abs(1.0 / (alpha + 1j * t)) ** 2)
    print(f"grid points        : {len(t)}")
    print(f"||E||^2 numerico   : {normE2:.6f}   (exacto 1/(2a) = {1/(2*alpha):.6f})")
    print()

    ns = [10, 20, 40, 80, 160, 320, 640, 1000, 1500]
    M_full, y, ks = assemble(alpha, t, w, max(ns))
    mu = np.array([mobius_function(k) for k in ks], dtype=float)

    print(f"{'n':>6} {'d_n^2 (optimo)':>16} {'s_n (mu(k))':>16} {'d_n^2*log n':>14} {'s_n*log n':>12}")
    for n in ns:
        j = n - 1                     # columnas k=2..n  ->  indices 0..n-2
        M = M_full[:, :j]
        # d_n^2 = residuo de minimos cuadrados
        c, *_ = np.linalg.lstsq(M, y, rcond=None)
        dn2 = float(np.real(np.vdot(y - M @ c, y - M @ c)))
        # s_n con coeficientes mu(k)
        r = M @ mu[:j] - y
        sn = float(np.real(np.vdot(r, r)))
        print(f"{n:>6} {dn2:>16.6e} {sn:>16.6e} {dn2*np.log(n):>14.6f} {sn*np.log(n):>12.4f}")

    print()
    print("Constante conjeturada de referencia: sum_rho 1/|rho|^2 = 2+gamma-log(4pi) ~ 0.0461")

if __name__ == "__main__":
    main()
