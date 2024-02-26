try:
    from oasys_srw.srwlib import *
    from oasys_srw.uti_plot import *
except:
    from srwlib import *
    from uti_plot import *

import numpy

import numpy

def _srw_stokes0_to_arrays(stk):
  Shape = (4,stk.mesh.ny,stk.mesh.nx,stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape,dtype=stk.arS.typecode)
  data0 = data[0]
  data1 = data[1]
  x = numpy.linspace(stk.mesh.xStart,stk.mesh.xFin,stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart,stk.mesh.yFin,stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart,stk.mesh.eFin,stk.mesh.ne)
  Z2 = numpy.zeros((e.size,x.size,y.size))
  POL_DEG = numpy.zeros((e.size,x.size,y.size))
  for ie in range(e.size):
      for ix in range(x.size):
          for iy in range(y.size):
            Z2[ie,ix,iy] = data0[iy,ix,ie]
            # this is shadow definition, that uses POL_DEG = |Ex|/(|Ex|+|Ey|)
            Ex = numpy.sqrt(numpy.abs(0.5*(data0[iy,ix,ie]+data1[iy,ix,ie])))
            Ey = numpy.sqrt(numpy.abs(0.5*(data0[iy,ix,ie]-data1[iy,ix,ie])))
            POL_DEG[ie,ix,iy] =  Ex / (Ex + Ey)
  return Z2,POL_DEG,e,x,y


def plot_and_write_files(wfr, root="id06_backpropagated"):
    # print(arI1.shape, plotMesh1x.shape, plotMesh1y.shape)

    stk = SRWLStokes()
    stk.mesh = mesh
    stk.allocate(mesh.ne, mesh.nx, mesh.ny)
    wfr.calc_stokes(stk)
    Z2, POL_DEG, e, x, y = _srw_stokes0_to_arrays(stk)
    print(Z2.shape, e.shape, x.shape, y.shape)
    from srxraylib.plot.gol import plot_image, plot

    plot_image(Z2[0], x, y)
    plot(x, Z2[0, :, y.size // 2],
         y, Z2[0, x.size // 2, :])

    z = Z2[0]

    filename = "%s_x.dat" % root
    f = open(filename, 'w')
    for i in range(y.size):
        f.write("%g  %g\n" % (x[i], z[i, y.size // 2]))
    f.close()
    print("File written to disk: %s" % filename)

    filename = "%s_y.dat" % root
    f = open(filename, 'w')
    for i in range(y.size):
        f.write("%g  %g\n" % (y[i], z[x.size // 2, i]))
    f.close()
    print("File written to disk: %s" % filename)



####################################################
# LIGHT SOURCE

part_beam = SRWLPartBeam()
part_beam.Iavg               = 0.2
part_beam.partStatMom1.x     = 0.0
part_beam.partStatMom1.y     = 0.0
part_beam.partStatMom1.z     = -1.071999
part_beam.partStatMom1.xp    = 0.0
part_beam.partStatMom1.yp    = 0.0
part_beam.partStatMom1.gamma = 11741.70710144324
part_beam.arStatMom2[0]      = 0.0
part_beam.arStatMom2[1]      = 0.0
part_beam.arStatMom2[2]      = 0.0
part_beam.arStatMom2[3]      = 0.0
part_beam.arStatMom2[4]      = 0.0
part_beam.arStatMom2[5]      = 0.0
part_beam.arStatMom2[10]     = 1.0000000000000001e-16

magnetic_fields = []
magnetic_fields.append(SRWLMagFldH(1, 'v',
                                   _B=0.7979326048108762,
                                   _ph=0.0,
                                   _s=-1,
                                   _a=1.0))
magnetic_structure = SRWLMagFldU(_arHarm=magnetic_fields, _per=0.018, _nPer=111.111)
magnetic_field_container = SRWLMagFldC(_arMagFld=[magnetic_structure],
                                       _arXc=array('d', [0.0]),
                                       _arYc=array('d', [0.0]),
                                       _arZc=array('d', [0.0]))

mesh = SRWLRadMesh(_eStart=10000.001989244456,
                   _eFin  =10000.001989244456,
                   _ne    =1,
                   _xStart=-0.00295,
                   _xFin  =0.00295,
                   _nx    =201,
                   _yStart=-0.00295,
                   _yFin  =0.00295,
                   _ny    =201,
                   _zStart=100.0)

stk = SRWLStokes()
stk.allocate(1,201,201)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh
wfr.partBeam = part_beam
wfr.unitElFld = 1

initial_mesh = deepcopy(wfr.mesh)
srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1,0.01,0.0,0.0,50000,1,0.0])

mesh0 = deepcopy(wfr.mesh)
arI = array('f', [0]*mesh0.nx*mesh0.ny)
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
arIx = array('f', [0]*mesh0.nx)
srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
arIy = array('f', [0]*mesh0.ny)
srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot_init(backend="Qt5Agg")
uti_plot2d1d (arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])

plot_and_write_files(wfr, root="id06_farfield")

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []

drift_before_oe_0 = SRWLOptD(-100.0)
pp_drift_before_oe_0 = [0,0,1.0,1,0,0.15,10.0,0.15,10.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_0)
srw_pp_array.append(pp_drift_before_oe_0)



####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)
srwl.PropagElecField(wfr, optBL)

mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny)
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0)
arI1x = array('f', [0]*mesh1.nx)
srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, 0)
arI1y = array('f', [0]*mesh1.ny)
srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI1, mesh1, <file_path>)
plotMesh1x = [1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx]
plotMesh1y = [1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny]
uti_plot2d1d(arI1, plotMesh1x, plotMesh1y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
uti_plot_show()

plot_and_write_files(wfr, root="id06_backpropagated")

