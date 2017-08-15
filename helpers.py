import numpy as np
import math
from astropy import wcs


def draw_ellipse(image, x, y, major, minor, pa):
    Y, X = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    X = X.ravel()
    Y = Y.ravel()

    # Add a half pixel to each of the major and minor axes
    # to cater to pixels being discrete in size
    #major += 0.5
    #minor += 0.5

    # Perform translation and rotation first
    Xdash = (X-x) * np.cos(np.radians(pa)) + (Y-y) * np.sin(np.radians(pa))
    Ydash = (Y-y) * np.cos(np.radians(pa)) - (X-x) * np.sin(np.radians(pa))

    ellipse = (Xdash**2 / major**2 + (Ydash)**2 / minor**2 <= 1)
    ellipse = ellipse.reshape((image.shape[0], image.shape[1]))

    image[ellipse] = 1


def elliptical_gaussian(X, Y, x, y, amp, major, minor, pa):
    """
    Generate a model 2d Gaussian with the given parameters.
    Evaluate this model at the given locations x,y.

    :param x,y: locations at which to calculate values
    :param xo,yo: position of Gaussian
    :param amp: amplitude of Gaussian
    :param major,minor: axes (sigmas)
    :param theta: position angle (degrees) CCW from x-axis
    :return: Gaussian function evaluated at x,y locations
    """
    try:
        sint, cost = math.sin(np.radians(pa)), math.cos(np.radians(pa))
    except ValueError as e:
        if 'math domain error' in e.args:
            sint, cost = np.nan, np.nan
        else:
            raise
    xxo = X - x
    yyo = Y - y
    exp = (xxo * cost + yyo * sint) ** 2 / (major) ** 2 + \
          (xxo * sint - yyo * cost) ** 2 / (minor) ** 2
    exp *= -1. / 2
    return amp * np.exp(exp)


def draw_gaussian(data, x, y, sx, sy, pa, peak):
    # Only apply Gaussian's out to this limit. Dramatically speeds up image creation.
    gaussian_limit = 6 * max(sx, sy)

    ysize, xsize = data.shape
    xmin = int(min(max(0, x - gaussian_limit), xsize))
    xmax = int(min(max(0, x + gaussian_limit), xsize))
    ymin = int(min(max(0, y - gaussian_limit), ysize))
    ymax = int(min(max(0, y + gaussian_limit), ysize))

    if xmin == xmax or ymin == ymax:
        return

    Y, X = np.mgrid[ymin:ymax, xmin:xmax]
    X = X.ravel()
    Y = Y.ravel()

    model = elliptical_gaussian(X, Y, x, y,
                                peak, sx,
                                sy, pa)

    # Blanking
    model = model.reshape((ymax-ymin, xmax-xmin))
    data[ymin:ymax, xmin:xmax] += model


def translate(ra, dec, r, theta):
    """
    Translate the point (ra,dec) a distance r (degrees)
    along angle theta (degrees)
    The translation is taken along an arc of a great circle.
    Return the (ra,dec) of the translated point.
    """
    factor = np.sin(np.radians(dec)) * np.cos(np.radians(r))
    factor += np.cos(np.radians(dec)) * np.sin(np.radians(r)) * np.cos(np.radians(theta))
    dec_out = np.degrees(np.arcsin(factor))

    y = np.sin(np.radians(theta)) * np.sin(np.radians(r)) * np.cos(np.radians(dec))
    x = np.cos(np.radians(r)) - np.sin(np.radians(dec)) * np.sin(np.radians(dec_out))
    ra_out = ra + np.degrees(np.arctan2(y, x))
    return ra_out, dec_out


# The following functions are explained
# at http://www.movable-type.co.uk/scripts/latlong.html
# phi ~ lat ~ Dec
# lambda ~ lon ~ RA
def gcd(ra1, dec1, ra2, dec2):
    """
    Great circle distance as calculated by the haversine formula
    ra/dec in degrees
    returns:
    sep in degrees
    """
    # TODO:  Vincenty formula
    # see - https://en.wikipedia.org/wiki/Great-circle_distance
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(np.radians(dlat) / 2) ** 2
    a += np.cos(np.radians(dec1)) * np.cos(np.radians(dec2)) * np.sin(np.radians(dlon) / 2) ** 2
    sep = np.degrees(2 * np.arcsin(min(1, np.sqrt(a))))
    return sep


class WCSHelper:
    def __init__(self, header):
        self.wcs = wcs.WCS(header)
        self.header = header

    def sky2pix(self, pos):
        """
        Take pos = (ra,dec) coords
        convert to pixel = (x,y) coords
        """
        x, y = pos
        # pixel = self.wcs.wcs_world2pix([(x, y, 0, 0)], 0)
        pixel = self.wcs.wcs_world2pix([(x, y)], 0)
        return [pixel[0][0], pixel[0][1]]

    def sky2pix_ellipse(self, pos, a, b, pa):
        """
        Convert an ellipse from sky to pixel corrds
        a/b vectors are calculated at an origin pos=(ra,dec)
        All input parameters are in degrees
        Output parameters are:
        x,y - the x,y pixels corresponding to the ra/dec position
        sx, sy - the major minor axes (FWHM) in pixels
        theta - the position angle in degrees

        :param pos: [ra,dec] of the ellipse center
        :param a: major axis
        :param b: minor axis
        :param pa: position angle
        :return: x, y, sx, sy, theta
        """
        ra, dec = pos
        x, y = self.sky2pix(pos)

        x_off, y_off = self.sky2pix(translate(ra, dec, a, pa))
        sx = np.hypot((x - x_off), (y - y_off))
        theta = np.arctan2((y_off - y), (x_off - x))

        x_off, y_off = self.sky2pix(translate(ra, dec, b, pa-90))
        sy = np.hypot((x - x_off), (y - y_off))
        theta2 = np.arctan2((y_off - y), (x_off - x)) - np.pi/2

        # The a/b vectors are perpendicular in sky space,
        # but not always in pixel space
        # so we have to account for this by calculating
        # the angle between the two vectors
        # and modifying the minor axis length
        defect = theta - theta2
        sy *= abs(np.cos(defect))

        return x, y, sx, sy, np.degrees(theta)

    def beamarea_pix(self):
        """
        Returns the area of the beam in pixels squared,
        taking into account the special 4 ln(2) factor.
        """
        beamsigma1 = self.header['BMAJ'] / self.wcs.wcs.cdelt[0]
        beamsigma2 = self.header['BMIN'] / self.wcs.wcs.cdelt[0]
        return (np.pi * beamsigma1 * beamsigma2) / (4 * np.log(2))