
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

    def __init__(nx, ny, distance, xangle, yangle, refpoint):

        self.nx = nx
        self.ny = ny
        self.distance = distance
        self.xangle = xangle
        self.yangle = yangle
        self.refpoint = refpoint


class Edge:

    def __init__(self, point1, point2):

        if point1.lon < point2.lon:
            wpoint = point1
            epoint = point2
        else:
            wpoint = point2
            epoint = point1

        self.lons = (wpoint.lon, epoint.lon)
        self.lats = (wpoint.lat, epoint.lat)

        cosd = cos(wpoint.lat)*cos(epoint.lat) + sin(wpoint.lat)*sin(epoint.lat)*cos(epoint.lon-wpoint.lon)

        self.cos = (cos(epoint.lat) - cos(wpoint.lat)*cosd)/(sin(wpoint.lat)*sin(acos(cosd)))


def compute_equidistant_point(lon1, lat1, lon2, lat2, d, right):
    """
    Determine a point that is distance d from point 1 and point 2
    If right is true, 
    the point will be to the right of the line between point 1 and point 2,
    otherwise, it will be to the left
    """

    if right:
        sign = 1
    else:
        sign = -1

    cosd3 = cos(lat1)*cos(lat2) + sin(lat1)*sin(lat2)*cos(lon1-lon2)
    sind3 = sin(acos(cosd3))

    costheta1p  = (cos(lat2) - cos(lat1)*cosd3)/(sin(lat1)*sind3)
    costheta1pp = (cos(d) - cos(d)*cosd3)/(sin(d)*sind3)

    theta1p  = acos(costheta1p)
    theta1pp = acos(costheta1pp)

    theta1 = acos(costheta1p) + sign*acos(costheta1pp)

    coslat3 = cos(lat1)*cos(d) + sin(lat1)*sin(d)*cos(theta1)
    lat3 = acos(coslat3)

    coslambda3 = (cos(d)-cos(lat1)*coslat3)/(sin(lat1)*sin(lat3))

    lon3 = min(lon1, lon2) + acos(coslambda3)

    return lon3, lat3


def generate_grid(nx, ny, distance, xangle, yangle, startpoint):

    points = [startpoint]
    for nn in range(1,ny):
        points += [points[-1].compute_point(distance, yangle)]

    for n in range(1,nx):
        points += [points[-ny].compute_point(distance, xangle)]
        for nn in range(1,ny):
            points += [points[-1].compute_point(distance, yangle)]


    for point in points:
        lat = point.get_deg_lat()
        lon = point.get_deg_lon()

        print("{:11.4f}{:11.4f}".format(lon, lat))



if __name__ == "__main__":

    nx = 7
    ny = 7

    distance = 2.5
    r  = 6371.0088

    anglex = 130
    angley = 90

    startlatp   = 71
    startlatpp  = 38
    startlatppp = 5

    startlonp   = 20
    startlonpp  = 40
    startlonppp = 31

    startpoint = Point.min_sec(startlatp, startlatpp, startlatppp, startlonp, startlonpp, startlonppp)

    xangle = Angle.degrees(anglex)
    yangle = Angle.degrees(anglex + angley)

    generate_grid(nx, ny, distance/r, xangle, yangle, startpoint)

