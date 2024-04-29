import vtk


def intp_k(x, points, const=2, clamp=True, tension=-0.3, bias=0, continuity=0):
    """
    Kochanek spline, wrapping `vtk.vtkKochanekSpline()
    <http://www.vtk.org/doc/nightly/html/classvtkKochanekSpline.html>`_
    args:
        x (list): t coordinate list (in range=[0,1])
        points (list) : coordinate x,y,z points
        const (int): constraint type for left and right of spline
        clamp (bool): flag on whether spline is clamped
        tension (float): tension parameter
        bias (float): bias parameter value
        continuity (int) : continuity parameter value
    returns:
        x : t coordinate list
        o : evaluated coordinate list
    """
    sc = vtk.vtkKochanekSpline()
    sc.SetLeftConstraint(const)
    sc.SetRightConstraint(const)
    sc.SetDefaultTension(tension)
    sc.SetDefaultBias(bias)
    sc.SetDefaultContinuity(continuity)
    if clamp:
        sc.ClampValueOn()
    for i in points:
        sc.AddPoint(i[0], i[1])

    o = [sc.Evaluate(i) for i in x]
    return x, o


def intp_c(x, points, const=2, clamp=True):
    """
    Cardinal spline, wrapping `vtk.vtkCardinalSpline()
    <http://www.vtk.org/doc/nightly/html/classvtkCardinalSpline.html>`_
    args:
        x (list): t coordinate list (in range=[0,1])
        points (list) : coordinate x,y,z points
        const (int): constraint type for left and right of spline
        clamp (bool): flag on whether spline is clamped
    """
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
        print(points)
        exit(e)


def intp_sc(x, points):
    """
    SCurve spline based interpolation
    args:
        x (list) : t coordinate list
        points (list) : xyz coordinate input points
    returns:
        x (relative coordinate point list)
        o (xyz coordinate points list, resplined)
    """
    sc = vtk.vtkSCurveSpline()
    for i in points:
        sc.AddPoint(i[0], i[1])
    o = [sc.Evaluate(i) for i in x]
    return x, o
