#! /usr/bin/env python3

import argparse
import csv
from pathlib import Path

from astropy.io import fits
import numpy as np

from helpers import draw_gaussian, WCSHelper


NVSS = open('/Users/torrance/ResearchB/Catalog/nvss.txt')
TGSS = open('/Users/torrance/ResearchB/Catalog/TGSSADR1_7sigma_catalog.tsv')
RATIO = 4.5


def main():
    parser = argparse.ArgumentParser(description="Blanks sources based on NVSS and TGSS catalogs")
    parser.add_argument('file', type=str)
    args = parser.parse_args()

    hdulist = fits.open(args.file)
    data = np.squeeze(hdulist[0].data)
    noise_mean, noise_std = noise(data)
    print("Noise: ", noise_mean, noise_std)

    # Extract beam size and coordinate information
    bmaj = hdulist[0].header['BMAJ']
    bmajpx = bmaj / abs(hdulist[0].header['CD1_1'])
    bmajsigma = bmajpx / 2.35482
    wcshelper = WCSHelper(hdulist[0].header)

    # Open output
    path = Path(args.file).with_suffix(".contour.fits")
    output_contour = path.open(mode='w')

    path = Path(args.file).with_suffix(".blanked.fits")
    output_blanked = path.open(mode='w')

    image = np.zeros(data.shape)

    tgss_sources = tgss_get_sources(TGSS)
    nvss_sources = nvss_get_sources(NVSS)
    sources = tgss_sources + nvss_sources
    for i, (ra, dec, a, b, wpa, peak) in enumerate(sources):
        print("Drawing source {}/{}".format(i, len(sources)), end="\r")
        # Convert WCS to pixel values
        x, y, sx, sy, pa = wcshelper.sky2pix_ellipse((ra, dec), a, b, wpa)

        # Convert FWHM values to sigmas
        sx = sx / 2.35482
        sy = sy / 2.35482

        # Convolve catalog source
        sx = np.sqrt(bmajsigma**2 + sx**2)
        sy = np.sqrt(bmajsigma**2 + sy**2)
        peak = (peak * 2 * np.pi) / np.sqrt((1/sx**2 + 1/bmajsigma**2) * (1/sy**2 + 1/bmajsigma**2))

        draw_gaussian(image, x, y, sx, sy, pa, peak)

    print("\nDone")

    hdu = fits.PrimaryHDU([[image]], header=hdulist[0].header.copy())
    hdu.writeto(output_contour, clobber=True)

    # 072-103: 0.5
    # 170 - 231: 0.01
    data[image > 3 * noise_std] = np.nan
    hdu = fits.PrimaryHDU([[data]], header=hdulist[0].header.copy())
    hdu.writeto(output_blanked, clobber=True)


def tgss_get_sources(tgss):
    sources = []
    reader = csv.reader(tgss, delimiter="\t")

    for i, row in enumerate(reader):
        # Skip header
        if i == 0:
            continue

        ra = float(row[1])
        dec = float(row[3])
        a = float(row[9]) / 3600  # Convert arcseconds -> degrees
        b = float(row[11]) / 3600  # Convert arcseconds -> degress
        wpa = float(row[13])
        peak = float(row[7]) / 1000  # mJy/beam

        # EoR0 field
        if not (ra < 35 or ra > 325) or not (-55 < dec < 0):
            continue

        # Spectral peak correction (ish!!)
        # peak = peak * RATIO
        sources.append([ra, dec, a, b, wpa, peak])

    return sources


def nvss_get_sources(nvss):
    sources = []
    for i, line in enumerate(nvss):
        # Skip headers
        if i < 17:
            continue

        parts = line.split()
        # Skip error lines and title fields which keep reappearing
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

        # Spectral correction (ish!!)
        peak = peak * RATIO
        sources.append([ra, dec, a, b, wpa, peak])

    return sources


def noise(image):
    image = image.ravel()
    for _ in range(3):
        mean = np.mean(image)
        std = np.std(image)
        image = image[image < mean + 3 * std]
        image = image[image > mean - 3 * std]

    return np.mean(image), np.std(image)

if __name__ == '__main__':
    main()
