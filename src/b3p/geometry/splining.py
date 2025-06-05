import vtk
import logging

logger = logging.getLogger(__name__)


def intp_c(x, points, const=2, clamp=True):
    try:
        sc = vtk.vtkCardinalSpline()
        sc.SetLeftConstraint(const)
        sc.SetRightConstraint(const)
        if clamp:
            sc.ClampValueOn()
        else:
            sc.ClampValueOff()
        for i in points:
            sc.AddPoint(i[0], i[1])

        o = [sc.Evaluate(i) for i in x]
        return x, o

    except Exception as e:
        logger.error(f"Error occurred while interpolating points: {points}")
        logger.exception(e)
