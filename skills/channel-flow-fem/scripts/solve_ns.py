"""Steady incompressible Navier-Stokes in a 2-D microchannel, Taylor-Hood P2/P1.

Physics: incompressible, Newtonian, steady, laminar, rigid walls. Picard iteration on the
convective term until the relative velocity change falls below tol.

Default fluid properties are those of a typical cell-culture medium at 37 C:
rho = 998 kg/m3, mu = 7e-4 Pa.s. Change RHO and MU for anything else.

The convective term is retained deliberately. "Laminar" does not mean "Stokes": at
Reynolds numbers of tens to a few hundred, cavities between ridges hold recirculating
vortices that a Stokes (reversible, inertia-free) model cannot produce at all.

2-D mid-plane approximation: flow rate is converted to mean velocity through the spanwise
width W_CHAN, so a result stands for the mid-plane of a channel whose width is much larger
than its height. Side walls are neglected -- the same parallel-plate idealisation behind the
closed-form wall shear tau = 6*mu*Q/(w*h^2) used to verify this solver.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import splu
from skfem import (MeshTri, Basis, FacetBasis, ElementVector, ElementTriP2, ElementTriP1,
                   BilinearForm, asm, condense, solve)
from skfem.helpers import grad, ddot, div

RHO, MU = 998.0, 7e-4
W_CHAN = 5e-3          # spanwise width, used only to convert flow rate -> mean velocity
H_CHAN = 500e-6        # inlet height


@BilinearForm
def vlap(u, v, w):
    return ddot(grad(u), grad(v))


@BilinearForm
def divu(u, q, w):
    return div(u) * q


@BilinearForm
def conv(u, v, w):
    return np.einsum("j...,ij...,i...->...", w["wind"], grad(u), v)


def inlet_profile(z, U_mean, H=H_CHAN):
    """Fully developed plane-Poiseuille profile, mean velocity U_mean."""
    return 1.5 * U_mean * (1.0 - ((z - 0.5 * H) / (0.5 * H)) ** 2)


def solve_case(mesh_path, Q_mL_min, H=H_CHAN, tol=1e-7, maxit=60, relax=1.0, verbose=True):
    m = MeshTri.load(mesh_path)
    ub = Basis(m, ElementVector(ElementTriP2()), intorder=4)
    pb = Basis(m, ElementTriP1(), intorder=4)
    Nu, Np = ub.N, pb.N

    Q = Q_mL_min * 1e-6 / 60.0
    U_mean = Q / (W_CHAN * H)

    A = asm(vlap, ub)
    B = asm(divu, ub, pb)
    Z = sp.csr_matrix((Np, Np))

    # Dirichlet: parabolic inlet, no-slip on both walls
    x = np.zeros(Nu + Np)
    din = ub.get_dofs("inlet")
    d1, d2 = din.all(["u^1"]), din.all(["u^2"])
    x[d1] = inlet_profile(ub.doflocs[1, d1], U_mean, H=H)
    x[d2] = 0.0
    D = np.unique(np.concatenate([din.all(), ub.get_dofs("top").all(), ub.get_dofs("bottom").all()]))

    uh = np.zeros(Nu)
    hist = []
    for it in range(maxit):
        wind = ub.interpolate(uh)
        C = asm(conv, ub, wind=wind)
        S = sp.bmat([[MU * A + RHO * C, -B.T], [B, Z]], format="csr")
        sol = solve(*condense(S, np.zeros(Nu + Np), x=x, D=D))
        u_new = sol[:Nu]
        den = np.linalg.norm(u_new) or 1.0
        err = np.linalg.norm(u_new - uh) / den
        uh = uh + relax * (u_new - uh) if relax != 1.0 else u_new
        hist.append(err)
        if verbose:
            print(f"    picard {it:2d}  rel.change {err:.3e}", flush=True)
        if err < tol:
            break
    else:
        raise RuntimeError(f"Picard did not converge: last {err:.2e}")

    ph = sol[Nu:]
    return dict(mesh=m, ub=ub, pb=pb, u=uh, p=ph, U_mean=U_mean, Q_mL_min=Q_mL_min,
                iters=it + 1, hist=hist, ndofs=Nu + Np)


def wall_shear(res, boundary="bottom"):
    """Signed tangential viscous traction on a wall, and the arclength/x of each point."""
    m, ub = res["mesh"], res["ub"]
    fb = FacetBasis(m, ub.elem, facets=m.boundaries[boundary])
    du = fb.interpolate(res["u"]).grad           # (i, j, nfacets, nqp) = du_i/dx_j
    n = fb.normals                               # (2, nfacets, nqp)
    # viscous traction t_i = mu (du_i/dx_j + du_j/dx_i) n_j
    t = MU * (np.einsum("ij...,j...->i...", du, n) + np.einsum("ji...,j...->i...", du, n))
    tn = np.einsum("i...,i...->...", t, n)
    # traction exerted BY the fluid ON the wall is -t (n is the outward normal of the fluid)
    tt = -(t - tn * n)                           # tangential part
    mag = np.sqrt(np.einsum("i...,i...->...", tt, tt))
    sign = np.sign(tt[0])                         # + when acting downstream
    xq = fb.global_coordinates().value[0]
    zq = fb.global_coordinates().value[1]
    return xq.ravel(), zq.ravel(), (sign * mag).ravel()


if __name__ == "__main__":
    import sys, json
    res = solve_case(sys.argv[1], float(sys.argv[2]))
    x, z, tau = wall_shear(res)
    print(json.dumps(dict(Q=res["Q_mL_min"], iters=res["iters"], ndofs=res["ndofs"],
                          tau_max_Pa=float(np.nanmax(tau)), tau_min_Pa=float(np.nanmin(tau)))))
