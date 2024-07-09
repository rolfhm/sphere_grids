from math import pi, cos, sin, acos
from itertools import pairwise
from sys import exit
from re import sub

r  = 6371.0088

class Point:

    def __init__(self, lat, lon):

        while lat < 0:
            lat += 2*pi
        while lat > 2*pi:
            lat -= 2*pi

        if lat > pi:
            lat = 2*pi - lat
            lon -= pi

        while lon < 0:
            lon += 2*pi
        while lon > 2*pi:
            lon -= 2*pi

        self.lat = lat
        self.lon = lon

    @classmethod
    def degrees(cls, lat, lon):
        """
        Initialise point based on latitude and longitude in degress
        """
        return cls(2*pi*(90 - lat)/360, 2*pi*lon/360)

    @classmethod
    def min_sec(cls, lat, latp, latpp, lon, lonp, lonpp):
        """
        Initialise point based on latitude and longitude in degrees, minutes, and seconds
        """
        return cls(2*pi*(90 - (lat + latp/60 + latpp/3600))/360, 2*pi*(lon + lonp/60 + lonpp/3600)/360)


    def compute_point(self, d3, theta1):
        """
        Compute new point (point 2) on a unit sphere (in radians)
        distance d3 from self (point 1)
        in the direction defined by theta1 between 0 and pi where 0.
        east is a logical, where True means the direction is east

        The point is computed using spherical trigonometry
        with an auxilliary point on the north pole (point 3)
        """

        if d3 == 0:
            return self

        if theta1.angle == 0:
            return Point(self.lat - d3, self.lon)

        if theta1.angle == pi:
            return Point(self.lat + d3, self.lon)

        d2 = self.lat

        cosd1 = cos(d3)*cos(d2) + sin(d3)*sin(d2)*cos(theta1.angle)
        d1 = acos(cosd1)

        costheta3 = (cos(d3) - cos(d1)*cos(d2))/(sin(d1)*sin(d2))
        theta3 = acos(costheta3)

        if theta1.east:
            return Point(d1, self.lon + theta3)
        else:
            return Point(d1, self.lon - theta3)


    def compute_distance(self, other):
        """
        Compute cosine of distance between self and other,
        using the north pole as a helping point
        """
        return cos(self.lat)*cos(other.lat) + sin(self.lat)*sin(other.lat)*cos(abs(self.lon-other.lon))


    def compute_cos(self, other):
        """
        Compute cosine of the angle between the north pole,
        self, and other
        cosd is the cosine of distance between self and other
        """

        if self.lon == other.lon:
            if self.lat < other.lat:
                return -1.0
            return 1.0

        cosd = self.compute_distance(other)
        return (cos(other.lat) - cos(self.lat)*cosd)/(sin(self.lat)*sin(acos(cosd)))


    def compute_angle(self, other):
        """
        Compute the angle between the north pole, self, and other
        """

        if self.lon == other.lon:
            if self.lat < other.lat:
                return Angle(0.0)
            return Angle(pi)

        cos = self.compute_cos(other)

        if self.lon < other.lon:
            return Angle(acos(cos))
        return Angle(2*pi - acos(cos))


    def lon_lat(self):
        """
        Print longitude and latitude of the gridpoint
        """
        return "{:11.4f}{:11.4f}".format(self.get_deg_lon(), self.get_deg_lat())


    def get_deg_lat(self):
        return 90 - 360*self.lat/(2*pi)

    def get_deg_lon(self):
        return 360*self.lon/(2*pi)

    def __eq__(self, other):
        return (self.lon == other.lon and self.lat == other.lat)

    def __repr__(self):
        return "({:8.6f}, {:8.6f})".format(self.lat, self.lon)


class Angle:
    """
    Store angle in normalized radians with logical indicating east or west.
    """

    def __init__(self, angle):
        """
        Compute east-west and store in objects.
        """

        while angle < 0:
            angle += 2*pi
        while angle > 2*pi:
            angle -= 2*pi

        if angle <= pi:
            self.angle = angle
            self.east = True
        else:
            self.angle = 2*pi - angle
            self.east = False


    @classmethod
    def degrees(cls, angle):
        """
        Initialise Angle object based on degrees
        """
        return cls(2*pi*angle/360)

    def get_abs(self):
        if self.east:
            return self.angle
        return 2*pi - self.angle

    def get_deg(self):
        return (360*self.get_abs())/(2*pi)

    def get_east(self):
        return self.east

    def add(self, other):
        """
        Add other to self and return new angle
        """
        return Angle(self.get_abs() + other.get_abs())

    def add_degree(self, angle):
        """
        Add other to self and return new angle
        """
        return self.add(Angle.degrees(angle))

    def __repr__(self):
        return str((self.angle, self.east))


class Grid:

    def __init__(self, nx, ny, dx, dy, xangle, yangle, refpoint):
        """
        Store the parameters of the grid and generate the gridpoints
        """
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.xangle = xangle
        self.yangle = yangle
        self.refpoint = refpoint

        self.generate_grid_points()


    def generate_grid_points(self):
        """
        Generate the points in the grid
        """
        points = [self.refpoint]
        for nn in range(1,self.ny):
            points += [points[-1].compute_point(self.dy, self.yangle)]

        for n in range(1,self.nx):
            points += [points[-self.ny].compute_point(self.dx, self.xangle)]
            for nn in range(1,self.ny):
                points += [points[-1].compute_point(self.dy, self.yangle)]

        self.points = points


    def filter(self, polygon, inside=True):
        """
        Filter out points in grid inside or outside the polygon
        depending on whether the points are inside or outside the polygon
        """

        filtered = []
        for point in self.points:
            if (polygon.inside(point) == inside):
                filtered += [point]

        self.points = filtered


    def lon_lat(self):
        """
        Return longitude and latitude of points
        """
        string = self.points[0].lon_lat()
        for point in self.points[1:]:
            string += '\n' + point.lon_lat()
        return string

        def __repr__(self):
            string = str(self.points[0])
            for point in self.points[1:]:
                string += '\n' + str(point)

        return(string)


class Edge:

    def __init__(self, pointx, pointy):
        """
        Initialise edge between pointx and pointy,
        store the westernmost as point1 and the easternmost point as point2,
        and store the cosine of the angle at the western point
        If the points are equally far west,
        store the northernmost point in point1
        """

        if pointx.lon == pointy.lon:
            self.aligned = True

            if pointx.lat < pointy.lat:
                self.point1 = pointx
                self.point2 = pointy
            else:
                self.point1 = pointy
                self.point2 = pointx

        else:
            self.aligned = False

            if pointx.lon < pointy.lon:
                self.point1 = pointx
                self.point2 = pointy
            else:
                self.point1 = pointy
                self.point2 = pointx

        self.cos1 = self.point1.compute_cos(self.point2)


    def lon1(self):
        return self.point1.lon

    def lon2(self):
        return self.point2.lon

    def lat1(self):
        return self.point1.lat

    def lat2(self):
        return self.point2.lat


    def point_crossing(self, point):
        """
        Check if a line between the north pole and point crosses the edge
        If point is on the edge, return True
        """
        if self.aligned:
            if (point.lon == self.lon1() and
                point.lat > self.lat1() and
                point.lat < self.lat2()):
                return True
            return False


        if (point.lon > self.lon1() and point.lon < self.lon2()):
            if (point.lat > self.lat1() and point.lat > self.lat2()):
                return True
            elif (point.lat > self.lat1() or point.lat > self.lat2()):
                if self.point1.compute_cos(point) < self.cos1:
                    return True
        return False


    def crossing(self, other):
        """
        Check if a self crosses with other
        First, check if longitudes and latitudes are overlapping,
        then check that the self.cos1 is between the angles
        from self.point1 to other.point1 and other.point2
        """
        if (self.point1 == other.point1 or
            self.point1 == other.point2 or
            self.point2 == other.point1 or
            self.point2 == other.point2):
            return False

        if (self.lon1() < other.lon2() and self.lon2() > other.lon1()):
            if (max(self.lat1(), self.lat2()) > min(other.lat1(), other.lat2()) and
                min(self.lat1(), self.lat2()) < max(other.lat1(), other.lat2())):

                    if (self.lon1() < other.lon1()):
                        angle1 = self.point1.compute_cos(other.point1)
                        angle2 = self.point1.compute_cos(other.point2)
                        cos1 = self.cos1
                    else:
                        angle1 = other.point1.compute_cos(self.point1)
                        angle2 = other.point1.compute_cos(self.point2)
                        cos1 = other.cos1

                    if (cos1 < max(angle1, angle2) and
                        cos1 > min(angle1, angle2)):
                        return True
        return False


    def lon_lat(self):
        """
        Print longitude and latitude of the gridpoints
        """
        return self.point1.lon_lat() + ',' + self.point2.lon_lat()


    def __repr__(self):
        return str(self.point1) + ', ' + str(self.point2)


class Polygon:

    def __init__(self, points):
        """
        Initialise the polygon based on corner points
        Store the edges in a list
        and the corners that need to be considered
        when determining if a point is inside the polygon or not
        """

        edges = []
        for point1, point2 in pairwise(points):
            edges += [Edge(point1, point2)]
        edges += [Edge(points[-1],points[0])]

        for i,edge in enumerate(edges[:-1]):
            for edgex in edges[i+1:]:
                if edge.crossing(edgex):
                    exit('crossing edges')

        self.edges = edges
        self.points = points

        crosspoints = []
        for i in range(len(points)):
            if ((points[i-1].lon < points[i].lon and points[i].lon < points[(i+1)%len(points)].lon) or
                (points[i-1].lon > points[i].lon and points[i].lon > points[(i+1)%len(points)].lon)):
                crosspoints += [points[i]]

        self.crosspoints = crosspoints


    @classmethod
    def min_sec(cls, latp, latpp, latppp, lonp, lonpp, lonppp):
        """
        Initialise polygon based on latitude and longitude in degrees, minutes, and seconds
        """
        points = []
        for ilatp, ilatpp, ilatppp, ilonp, ilonpp, ilonppp in zip(latp, latpp, latppp, lonp, lonpp, lonppp):
            points += [Point.min_sec(ilatp, ilatpp, ilatppp, ilonp, ilonpp, ilonppp)]

        return cls(points)


    def inside(self, point):
        """
        Determine if point is inside or outside self by checking
        if a line from the north pole crosses edges and corners
        an odd or even number of times
        """
        crossings = 0
        for edge in self.edges:
            if edge.point_crossing(point):
                crossings += 1

        for crosspoint in self.crosspoints:
            if (point.lon == crosspoint.lon and point.lat >= crosspoint.lat):
                crossings += 1

        if crossings%2 == 0:
            return False
        return True

    def lon_lat(self):
        """
        Print longitude and latitude of the points
        """
        string = self.points[0].lon_lat()
        for point in self.points[1:]:
            string += '\n' + point.lon_lat()
        return string

    def __repr__(self):
        string = str(self.edges[0])
        for edge in self.edges[1:]:
            string += '\n' + str(edge)

        return(string)

    def Diana_print(self):
        """
        Print points in a format easy to paste into Diana .kml files
        """
        string = ""
        for point in self.points:
            string += sub("\s+", ",", point.lon_lat().strip())
            string += "\n"
        return string

