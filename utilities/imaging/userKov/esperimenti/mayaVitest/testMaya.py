def dataGen():
    ''' generate the data in test.vtk file.'''
    from Numeric import *
    import scipy
    x = (arange(50.0)-25)/2.0
    y = (arange(50.0)-25)/2.0
    r = sqrt(x[:,NewAxis]**2+y**2)
    z = 5.0*scipy.special.j0(r) # Bessel function of order 0
    # now dump the data to a VTK file.
    import pyvtk
    # Flatten the 2D array data as per VTK's requirements.
    z1 = reshape(transpose(z), (-1,))
    point_data = pyvtk.PointData(pyvtk.Scalars(z1)) 
    grid = pyvtk.StructuredPoints((50,50, 1), (-12.5, -12.5, 0), (0.5, 0.5, 1)) 
    data = pyvtk.VtkData(grid, point_data) 
    data.tofile('test.vtk')

def dataVis():
    import mayavi
    v = mayavi.mayavi() # create a MayaVi window.
    d = v.open_vtk('test.vtk', config=0) # open the data file.
    # The config option turns on/off showing a GUI control for the data/filter/module.
    # load the filters.
    f = v.load_filter('WarpScalar', config=0)
    n = v.load_filter('PolyDataNormals', 0)
    n.fil.SetFeatureAngle (45) # configure the normals.
    # Load the necessary modules.
    m = v.load_module('SurfaceMap', 0)
    a = v.load_module('Axes', 0)
    a.axes.SetCornerOffset(0.0) # configure the axes module.
    o = v.load_module('Outline', 0)
    v.Render() # Re-render the scene.

if __name__=="__main__":
    #dataGen()
    dataVis()
