
from math import pi, cos, sin, acos

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


    def compute_angle(self, other):
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


    def get_deg_lat(self):
        return 90 - 360*self.lat/(2*pi)

    def get_deg_lon(self):
        return 360*self.lon/(2*pi)

    def __repr__(self):
        return str((self.lat, self.lon))


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

    def get_angle_deg(self):
        return 360*self.lon/(2*pi)

    def get_east(self):
        return self.east

    def __repr__(self):
        return str((self.angle, self.east))


class Grid:

    def __init__(self, nx, ny, distance, xangle, yangle, refpoint):
        """
        Store the parameters of the grid and generate the gridpoints
        """
        self.nx = nx
        self.ny = ny
        self.distance = distance
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
            points += [points[-1].compute_point(self.distance, self.yangle)]

        for n in range(1,self.nx):
            points += [points[-self.ny].compute_point(self.distance, self.xangle)]
            for nn in range(1,self.ny):
                points += [points[-1].compute_point(self.distance, self.yangle)]

        self.points = points


    def print_lat_lon(self):
        """
        Print latitude and longitude of the gridpoints
        """
        for point in self.points:
            lat = point.get_deg_lat()
            lon = point.get_deg_lon()

            print("{:11.4f}{:11.4f}".format(lat, lon))




class Edge:

    def __init__(self, pointx, pointy):
        """
        Initialise edge between pointx and pointy, 
        store the westernmost as point_1 and the easternmost point as point_2,
        and store the cosine of the angle at the western point
        If the points are equally far west,
        store the northernmost point in point_1
        """

        if pointx.lon == pointy.lon:
            self.aligned = True

            if pointx.lat < pointy.lat:
                self.point_1 = pointx
                self.point_2 = pointy
            else:
                self.point_1 = pointy
                self.point_2 = pointx

        else:
            self.aligned = False

            if pointx.lon < pointy.lon:
                self.point_1 = pointx
                self.point_2 = pointy
            else:
                self.point_1 = pointy
                self.point_2 = pointx

        self.cos_1 = self.point_1.compute_angle(self.point_2)


    def point_crossing(self, point):
        """
        Check if a line between the north pole and point crosses the edge
        If point is on the edge, return True
        """
        if self.aligned:
            if (point.lon == self.point_1.lon and
                point.lat > self.point_1.lat and 
                point.lat < self.point_2.lat):
                return True
            return False
            

        if (point.lon > self.point_1.lon or point.lon < self.point_2.lon):
            if (point.lat > self.point_1.lat and point.lat > self.point_2.lat):
                return True
            elif (point.lat > self.point_1.lat or point.lat > self.point_2.lat):
                if self.point_1.compute_angle(point) < self.cos_1:
                    print()
                    print(self.point_1.lon)
                    print(point.lon)
                    print(self.point_1.compute_angle(point))
                    print(self.cos_1)
                    print()
                    return True
        return False


if __name__ == "__main__":

    nx = 3
    ny = 3

    distance = 250
    r  = 6371.0088

    anglex = 135
    angley = 90

    startlatp   = 71
    startlatpp  = 38
    startlatppp = 5

#    startlonp   = 20
#    startlonpp  = 40
#    startlonppp = 31

    startlonp   = 0
    startlonpp  = 0
    startlonppp = 0

    startpoint = Point.min_sec(startlatp, startlatpp, startlatppp, startlonp, startlonpp, startlonppp)

    xangle = Angle.degrees(anglex)
    yangle = Angle.degrees(anglex + angley)

    grid = Grid(nx, ny, distance/r, xangle, yangle, startpoint)

#    grid.print_lat_lon()

    f_point = grid.points[0]
    l_point = grid.points[-1]

    corner1 = grid.points[6]
    corner2 = grid.points[-3]

    center = grid.points[4]

    edge = Edge(f_point, corner1)

    print()
    for point in grid.points:
        if point.lon > pi:
            print("{:11.4f}{:11.4f}".format(point.lat, point.lon - 2*pi))
        else:
            print("{:11.4f}{:11.4f}".format(point.lat, point.lon))
        
    print()
    print(edge.point_1, edge.point_1.get_deg_lat(), edge.point_1.get_deg_lon())
    print(edge.point_2, edge.point_2.get_deg_lat(), edge.point_2.get_deg_lon())
    print(edge.cos_1)
    print(acos(edge.cos_1), 360*acos(edge.cos_1)/(2*pi))
    print()
    print(edge.point_crossing(l_point))
    print(edge.point_crossing(corner2))
    print(edge.point_crossing(center))
    print()



