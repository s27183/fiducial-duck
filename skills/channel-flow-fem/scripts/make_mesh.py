"""Ridged microchannel geometry for gmsh, as a streamwise-vertical (x-z) 2-D section.

SI units (metres). build() writes a .msh file with physical groups "fluid", "inlet",
"outlet", "bottom" and "top"; n_ridge=0 gives a plain parallel-plate channel.

Defaults
--------
channel height                  : 500 um
ridge height x streamwise width : 100 um x 100 um, on the bottom wall
inter-ridge gap                 : 400 um   -> pitch 500 um
mesh size                       : lc_bulk in the core, lc_wall at walls and ridge faces

Modelling choices the caller should state when reporting results
----------------------------------------------------------------
* the number of ridges, and which cavity along the channel is analysed
* the straight inlet and outlet lengths; solve_ns imposes a fully developed inlet profile
  rather than resolving the entrance length
* the 2-D section stands for the mid-plane of a channel much wider than it is tall
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gmsh, numpy as np, sys

H       = 500e-6      # channel height (cavity region)
HR      = 100e-6      # ridge height
WR      = 100e-6      # ridge streamwise width
GAP     = 400e-6      # inter-ridge gap
N_RIDGE = 10
L_IN    = 1.5e-3
L_OUT   = 2.0e-3

PITCH = WR + GAP
L = L_IN + N_RIDGE * PITCH + L_OUT

lc_bulk = 25e-6
lc_wall = 6e-6


def build(path="ridged_channel.msh", n_ridge=N_RIDGE, lc_bulk=lc_bulk, lc_wall=lc_wall,
          L_IN=L_IN, L_OUT=L_OUT, H=H):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ridged")
    occ = gmsh.model.occ

    total_L = L_IN + n_ridge * PITCH + L_OUT
    chan = occ.addRectangle(0, 0, 0, total_L, H)
    ridges = []
    for i in range(n_ridge):
        x0 = L_IN + i * PITCH
        ridges.append((2, occ.addRectangle(x0, 0, 0, WR, HR)))
    if ridges:
        out, _ = occ.cut([(2, chan)], ridges)
    else:
        out = [(2, chan)]
    occ.synchronize()

    surf = [s[1] for s in out]

    # classify boundary curves
    inlet, outlet, bottom, top = [], [], [], []
    for dim, tag in gmsh.model.getEntities(1):
        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        x, y = com[0], com[1]
        if abs(x) < 1e-9:
            inlet.append(tag)
        elif abs(x - total_L) < 1e-9:
            outlet.append(tag)
        elif abs(y - H) < 1e-9:
            top.append(tag)
        else:
            bottom.append(tag)

    gmsh.model.addPhysicalGroup(2, surf, 1);    gmsh.model.setPhysicalName(2, 1, "fluid")
    gmsh.model.addPhysicalGroup(1, inlet, 1);   gmsh.model.setPhysicalName(1, 1, "inlet")
    gmsh.model.addPhysicalGroup(1, outlet, 2);  gmsh.model.setPhysicalName(1, 2, "outlet")
    gmsh.model.addPhysicalGroup(1, bottom, 3);  gmsh.model.setPhysicalName(1, 3, "bottom")
    gmsh.model.addPhysicalGroup(1, top, 4);     gmsh.model.setPhysicalName(1, 4, "top")

    # refine towards the ridged (bottom) wall where the shear peaks are reported
    df = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(df, "CurvesList", bottom)
    tf = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(tf, "InField", df)
    gmsh.model.mesh.field.setNumber(tf, "SizeMin", lc_wall)
    gmsh.model.mesh.field.setNumber(tf, "SizeMax", lc_bulk)
    gmsh.model.mesh.field.setNumber(tf, "DistMin", 20e-6)
    gmsh.model.mesh.field.setNumber(tf, "DistMax", 200e-6)
    gmsh.model.mesh.field.setAsBackgroundMesh(tf)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.mesh.generate(2)
    gmsh.write(path)
    nodes = gmsh.model.mesh.getNodes()[0].size
    _, elts = gmsh.model.mesh.getElementsByType(2)[0], gmsh.model.mesh.getElementsByType(2)[1]
    ncell = elts.size // 3
    gmsh.finalize()
    return dict(path=path, L=total_L, H=H, n_ridge=n_ridge, nodes=int(nodes), tris=int(ncell),
                lc_bulk=lc_bulk, lc_wall=lc_wall, pitch=PITCH)


if __name__ == "__main__":
    import json
    print(json.dumps(build(), indent=1))
