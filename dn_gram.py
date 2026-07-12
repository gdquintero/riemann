"""
Metodo acelerado para d_n(alpha)^2 via matriz de Gram y transformada de Fourier.

Productos internos en H^2(C_alpha):  <F,H> = (1/2pi) int_R F(a+it) conj(H(a+it)) dt.
Con G_k(a+it) = a_k(t) psi(t),  a_k(t) = k^{-a} e^{-it log k} - 1/k,  psi = zeta(a+it)/(a+it).

Peso:  rho(t) = |psi(t)|^2 = |zeta(a+it)|^2 / (a^2 + t^2)   (par y real)
Transformada:  W(x) = (1/2pi) int rho(t) e^{itx} dt   (real, par)

  A_{jk} = (jk)^{-a} W(log k - log j) - (j^{-a}/k) W(log j)
           - (k^{-a}/j) W(log k) + (1/(jk)) W(0)

Vector b_k = <E,G_k>, E=1/(a+it):
  g(t) = conj(psi(t))/(a+it) = conj(zeta(a+it))/(a^2+t^2)   (hermitiana: g(-t)=conj g(t))
  V(x) = (1/2pi) int g(t) e^{itx} dt   (real)
  b_k  = k^{-a} V(log k) - (1/k) V(0)

||E||^2 = 1/(2a).   d_n^2 = ||E||^2 - b^T A^{-1} b.

W, V se obtienen por FFT de rho, g sobre una malla uniforme en t.
"""
import numpy as np
from scipy.special import zeta

def transforms(alpha, dt=0.02, T=20000.0):
    """Devuelve (xgrid, W, V, W0, V0) con W,V muestreadas en xgrid = 2pi p/(M dt)."""
    t = np.arange(0.0, T, dt)
    M = len(t)
    s = alpha + 1j * t
    Z = zeta(s)
    denom = alpha**2 + t**2
    rho = (np.abs(Z) ** 2) / denom          # |psi|^2, real par
    g = np.conjugate(Z) / denom             # g(t) = conj(zeta)/(a^2+t^2)
    # pesos trapezoidales (extremo t=0 a la mitad; cola en T despreciable)
    cw = np.ones(M); cw[0] = 0.5
    # R(x_p) = dt * sum c_m f_m e^{+i 2pi m p/M} = dt * M * ifft(c f)[p]
    # F(x_p) = (1/pi) Re R.  (ifft lleva el signo +i correcto para f complejo)
    Wfull = (dt / np.pi) * M * np.real(np.fft.ifft(cw * rho))
    Vfull = (dt / np.pi) * M * np.real(np.fft.ifft(cw * g))
    P = M // 2
    xgrid = 2.0 * np.pi * np.arange(P) / (M * dt)
    return xgrid, Wfull[:P], Vfull[:P], float(Wfull[0]), float(Vfull[0])

def make_interp(xgrid, F):
    return lambda x: np.interp(np.abs(x), xgrid, F)   # W,V pares -> |x|

def gram(alpha, ns, dt=0.02, T=20000.0, eig_tol=1e-12):
    xgrid, Wg, Vg, W0, V0 = transforms(alpha, dt, T)
    Wf = make_interp(xgrid, Wg)
    Vf = make_interp(xgrid, Vg)
    normE2 = 1.0 / (2.0 * alpha)

    nmax = max(ns)
    ks = np.arange(2, nmax + 1)
    L = np.log(ks)
    Kam = ks.astype(float) ** (-alpha)
    Kinv = 1.0 / ks
    WL = Wf(L)                                  # W(log k)
    VL = Vf(L)                                  # V(log k)

    # matriz completa (para nmax moderado); si es grande, ojo con memoria
    Dl = np.abs(L[None, :] - L[:, None])        # |log k - log j|
    Wmat = Wf(Dl)
    A_full = (np.outer(Kam, Kam) * Wmat
              - np.outer(Kam * WL, Kinv)
              - np.outer(Kinv, Kam * WL)
              + np.outer(Kinv, Kinv) * W0)
    b_full = Kam * VL - Kinv * V0

    out = {}
    for n in ns:
        j = n - 1
        A = A_full[:j, :j]
        b = b_full[:j]
        # d_n^2 = ||E||^2 - b^T A^+ b  (A simetrica PSD, truncar autovalores chicos)
        ev, U = np.linalg.eigh(A)
        thr = eig_tol * ev.max()
        inv = np.where(ev > thr, 1.0 / ev, 0.0)
        proj = float((U.T @ b) ** 2 @ inv)
        out[n] = normE2 - proj
    return out, normE2

# ----------------- validacion contra metodo denso -----------------
def _dense_reference(alpha, nval, T_fine=60.0, h_fine=0.02, T_max=2.0e5, n_geo=4000):
    from proto_dn import build_grid, assemble
    t, w = build_grid(T_fine, h_fine, T_max, n_geo)
    M, y, ks = assemble(alpha, t, w, nval)
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    r = M @ c - y
    return float(np.real(np.vdot(r, r)))

def validate():
    alpha = 0.5
    ns = [40, 80, 160, 320]
    out, _ = gram(alpha, ns, dt=0.02, T=20000.0)
    print("Validacion Gram(FFT) vs denso (lstsq):")
    print(f"{'n':>5} {'Gram':>13} {'denso':>13} {'dif rel':>10}")
    for n in ns:
        dref = _dense_reference(alpha, n)
        rel = abs(out[n] - dref) / dref
        print(f"{n:>5} {out[n]:>13.6e} {dref:>13.6e} {rel:>10.2e}")

if __name__ == "__main__":
    validate()
