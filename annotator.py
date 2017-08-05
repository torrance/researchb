#! /usr/bin/env python3

import csv
import math
import numpy as np
from astropy.io import fits
from astropy import wcs


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

    Y, X = np.mgrid[ymin:ymax, xmin:xmax]
    X = X.ravel()
    Y = Y.ravel()

    model = elliptical_gaussian(X, Y, x, y,
                                peak, sx,
                                sy, pa)

    # Blanking
    model[model > 0.0001] = np.nan
    model[~np.isnan(model)] = 0
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
        pixel = self.wcs.wcs_world2pix([(x, y, 0, 0)], 0)
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


# Grab the header info
#hdulist = fits.open('../Images/TGSSADR_R01D19_5x5.MOSAIC.FITS')
hdulist = fits.open('../Images/EoR0-full-integration-for-Tara.fits')
wcshelper = WCSHelper(hdulist[0].header)
data = np.squeeze(hdulist[0].data).copy()

# print(repr(hdulist[0].header))

# Convolve the image with the larger beam
bmaj = 148.3 / 3600  # Median value from Aegean
# bmaj = 141 / 3600  # 25th percentile
# bmaj = 159 / 3600  # 75th percentile
# bmaj = 133.9 / 3600  # 10th percentile
# bmaj = 125 / 3600
print("{} versus header {}".format(bmaj, hdulist[0].header['BMAJ']))
# bmaj = hdulist[0].header['BMAJ']
bmajpx = bmaj / abs(hdulist[0].header['CDELT1'])
bmajsigma = bmajpx / 2.35482

sources = []

annf = open('../Output/tgss.ann', mode='w')
with open('../Catalog/TGSSADR1_7sigma_catalog.tsv') as f:
    reader = csv.reader(f, delimiter="\t")

    for i, row in enumerate(reader):
        # Skip header
        if i == 0:
            continue

        ra = float(row[1])
        dec = float(row[3])
        emaj = float(row[9]) / 3600  # Convert arcseconds -> degrees
        emin = float(row[11]) / 3600  # Convert arcseconds -> degress
        wpa = float(row[13])
        peak = float(row[7]) / 1000  # mJy/beam

        # if not (ra < 2.5 or ra > 357.5) or not (-24 < dec < -18):
        if not (ra < 35 or ra > 325) or not (-55 < dec < 0):
            continue

        if row[16] == 'S':
            print("COLOR GREEN", file=annf)
            print("ELLIPSE W {ra} {dec} {maj} {min} {pa}".format(
                ra=ra,
                dec=dec,
                min=emaj,
                maj=emin,
                pa=-wpa
            ), file=annf)
        if row[16] == 'M':
            print("COLOR ORANGE", file=annf)
            print("ELLIPSE W {ra} {dec} {maj} {min} {pa}".format(
                ra=ra,
                dec=dec,
                min=emaj,
                maj=emin,
                pa=-wpa
            ), file=annf)
        if row[16] == 'C':
            print("COLOR PINK", file=annf)
            print("ELLIPSE W {ra} {dec} {maj} {min} {pa}".format(
                ra=ra,
                dec=dec,
                min=emaj,
                maj=emin,
                pa=-wpa
            ), file=annf)

        # Convert WCS to pixel values
        x, y, sx, sy, pa = wcshelper.sky2pix_ellipse((ra, dec), emaj, emin, wpa)
        # Convert FWHM values to sigmas
        sx = sx / 2.35482
        sy = sy / 2.35482
        # The following transformations are the result a convolution
        # with the beam
        sx = np.sqrt(bmajsigma**2 + sx**2)
        sy = np.sqrt(bmajsigma**2 + sy**2)
        peak = (peak * 2 * np.pi) / np.sqrt((1/sx**2 + 1/bmajsigma**2) * (1/sy**2 + 1/bmajsigma**2))

        # Spectral peak correction (ish!!)
        peak = peak * 0.0906198
        sources.append([x, y, sx, sy, pa, peak])

annf.close()


annf = open('../Output/nvss.ann', 'w')
print("COLOR PINK", file=annf)
with open('../Catalog/nvss.txt') as f:
    for i, line in enumerate(f):
        # Skip header
        if i < 17:
            continue

        parts = line.split()

        if len(parts) < 14:
            continue

        # RA and Dec into degrees
        ra = float(parts[0]) * (360 / 24) + float(parts[1]) * (360 / (24 * 60)) + float(parts[2]) * (360 / (24 * 3600))
        dec = abs(float(parts[3])) + float(parts[4]) / 60 + float(parts[5]) / 3600
        if float(parts[3]) < 0:
            dec = -dec

        peak = float(parts[7]) / 1000
        a = float(parts[8]) / 3600
        b = float(parts[9]) / 3600
        wpa = float(parts[10])

        print("ELLIPSE W {ra} {dec} {maj} {min} {pa}".format(
            ra=ra,
            dec=dec,
            maj=a,
            min=b,
            pa=wpa
        ), file=annf)

        # Convert WCS to pixel values
        x, y, sx, sy, pa = wcshelper.sky2pix_ellipse((ra, dec), a, b, wpa)
        # Convert FWHM values to sigmas
        sx = sx / 2.35482
        sy = sy / 2.35482
        # The following transformations are the result a convolution
        # with the beam
        sx = np.sqrt(bmajsigma**2 + sx**2)
        sy = np.sqrt(bmajsigma**2 + sy**2)
        peak = (peak * 2 * np.pi) / np.sqrt((1/sx**2 + 1/bmajsigma**2) * (1/sy**2 + 1/bmajsigma**2))

        # Spectral correction (ish!!)
        peak = peak * 0.04119
        sources.append([x, y, sx, sy, pa, peak])

annf.close()

data0 = data.copy()
data0[:] = 0

for (x, y, sx, sy, pa, peak) in sources:
    draw_gaussian(data0, x, y, sx, sy, pa, peak)

# Blank out edges
# data0[:, :700] = np.nan
# data0[:, 4350:] = np.nan
# data0[:1072, :] = np.nan
# data0[4450:, :] = np.nan

# multipliers = []
# for i, (x, y, _, _, _, _) in enumerate(sources):
#     x = int(x)
#     y = int(y)
#     sources[i].append(1)  # Default value
#     if (620 < x < 4420) and (1340 < y < 4560):
#     #if (1 < x < data.shape[1]-2) and (1 < y < data.shape[0]-2):
#         if np.mean(data[y-1:y+1, x-1:x+1]) > 10:
#             multiplier = data[y-1:y+1, x-1:x+1] / data0[y-1:y+1, x-1:x+1]
#             multiplier = np.mean(multiplier)
#             sources[i][6] = multiplier
#             multipliers.append(multiplier)

# print("Multipliers:")
# print(np.mean(multipliers), np.std(multipliers))

# data0 = data.copy()
# data0[:] = 0
# for (x, y, sx, sy, pa, peak, multiplier) in sources:
#     draw_gaussian(data0, x, y, sx, sy, pa, peak*multiplier)

result = data.copy() - data0

# fig = FITSFigure(hdulist[0])
# fig.recenter(0, -25, radius=1.5)
# fig.show_colorscale(cmap='viridis')
# fig.save('original.pdf', dpi=2000)

# hdulist[0].data = [[data0]]
# hdulist.writeto('catalog-raw.fits', clobber=True)

# fig = FITSFigure(hdulist[0])
# fig.recenter(0, -25, radius=1)
# fig.show_colorscale(cmap='viridis')
# fig.save('catalog-raw.pdf', dpi=2000)

# hdulist[0].data = [[data0]]
# hdulist.writeto('catalog-convolved.fits', clobber=True)

# hdulist[0].data = [[data0]]
# hdulist.writeto('catalog-TGSS.fits', clobber=True)

# fig = FITSFigure(hdulist[0])
# fig.recenter(0, -25, radius=1.5)
# fig.show_colorscale(cmap='viridis')
# fig.save('catalog-convolved.pdf', dpi=2000)

# hdulist[0].data = [[result]]
# hdulist.writeto('subtracted.fits', clobber=True)

hdulist[0].data = [[result]]
hdulist.writeto('blanked.fits', clobber=True)

# fig = FITSFigure(hdulist[0])
# fig.recenter(0, -25, radius=1.5)
# fig.show_colorscale(cmap='viridis')
# fig.save('subtracted.pdf', dpi=2000)


