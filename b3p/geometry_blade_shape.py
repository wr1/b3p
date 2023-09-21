import vtk
from b3p import geometry_section
import numpy
from b3p import geom_utils


class blade_shape:
    def __init__(
        self,
        sections=None,
        section_resolution=200,
        web_resolution=20,
        added_datums=None,
    ):
        if sections is None:
            sections = []
        if added_datums is None:
            added_datums = {}
        self.set_sections(sections)
        self.set_section_resolution(section_resolution)
        self.webs = []
        self.web_resolution = web_resolution
        self.added_datums = added_datums

    def set_web(self, web):
        self.webs.append(web)

    def set_sections(self, sections):
        self.sections = sections

    def set_section_resolution(self, n_points):
        "set the number of points used to sample the airfoil section"
        self.n_points = n_points

    def build_interpolated_sections(self, radii, interpolation_type=1):
        """
        Create interpolated sections

        args:
            interpolation_type (int) : 1==cardinal spline, 2==linear, 3=kochanekspline
        """
        r = [i.r for i in self.sections]
        # get the resplined sections
        pnts = [i.respline(self.n_points)[0] for i in self.sections]

        # define n_points splines along the blade
        nxyz = []
        for i in zip(*pnts):
            if interpolation_type == 1:
                nx = geom_utils.spline_interp(r, list(zip(*i))[0], radii)
                ny = geom_utils.spline_interp(r, list(zip(*i))[1], radii)
                nz = geom_utils.spline_interp(r, list(zip(*i))[2], radii)
            elif interpolation_type == 2:
                nx = numpy.interp(radii, r, list(zip(*i))[0])
                ny = numpy.interp(radii, r, list(zip(*i))[1])
                nz = numpy.interp(radii, r, list(zip(*i))[2])
            elif interpolation_type == 3:
                nx = geom_utils.spline_interp_k(r, list(zip(*i))[0], radii)
                ny = geom_utils.spline_interp_k(r, list(zip(*i))[1], radii)
                nz = geom_utils.spline_interp_k(r, list(zip(*i))[2], radii)

            nxyz.append(list(zip(nx, ny, nz)))

        self.interp_sections = [
            geometry_section.section(i[0], i[0] / max(radii), i[1])
            for i in zip(radii, zip(*nxyz))
        ]

    def mesh(self, n_points=100, close=True, panel_mesh_scale=None):
        """
        Mesh generation for blade

        args:
            n_points (int) : number of points around cross section

            close (bool) : close two ends of mesh (i.e. TE)
        """
        if panel_mesh_scale is None:
            panel_mesh_scale = []
        self.poly = vtk.vtkPolyData()
        points = vtk.vtkPoints()

        # this dict is used to create a number of vtkFloatArrays that have named
        # variables as defined on the section, these can be used to define ply
        # positions
        added_arrays = {}
        # loop over the interpolated sections
        for i in self.interp_sections:
            # respline the section with the standard number of points
            pnts, data = i.respline(
                n_points, self.webs, self.added_datums, panel_mesh_scale
            )
            for j in pnts:
                points.InsertNextPoint(j[0], j[1], j[2])

            for j in data:
                if j not in added_arrays:
                    added_arrays[j] = vtk.vtkFloatArray()
                    added_arrays[j].SetName(j)
                for k in data[j]:
                    added_arrays[j].InsertNextValue(k)

        self.poly.SetPoints(points)

        np = self.poly.GetNumberOfPoints()

        quads = vtk.vtkCellArray()
        for i in range(1, int(np / n_points)):
            np = range((i - 1) * n_points, i * n_points)  # previous row point ids
            nc = range(i * n_points, (i + 1) * n_points)  # current row point ids
            for j in range(n_points - (0 if close else 1)):
                quads.InsertNextCell(4)
                quads.InsertCellPoint(np[j])
                quads.InsertCellPoint(nc[j])
                quads.InsertCellPoint(nc[(j + 1) % n_points])
                quads.InsertCellPoint(np[(j + 1) % n_points])

        self.poly.SetPolys(quads)
        for value in added_arrays.values():
            self.poly.GetPointData().AddArray(value)
        for i in self.webs:
            i.mesh(self.poly, self.web_resolution)

    def write_mesh(self, filename):
        """
        Export mesh in vtp format

        args:
            filename (str): vtp output filename
        """
        try:
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(filename)
            writer.SetInputData(self.poly)
            writer.Write()
        except Exception:
            print("no valid mesh available")

        for i in self.webs:
            i.write_mesh(f'{i.name.replace(".txt", ".vtp")}')
