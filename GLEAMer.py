#! /usr/bin/env python3

import csv
import pathlib
import os

from astropy.io import fits
import numpy as np

from helpers import WCSHelper, draw_ellipse


CATALOG = '/Users/torrance/Google Drive/Research B/Output/blanked_comp.csv'
FOLDER = '/Users/torrance/Google Drive/Research B/Stamps'

# islands = { island_id: { sources: [], stats: {} } }
islands = {}
for i, row in enumerate(csv.reader(open(CATALOG))):
    # Skip headers
    if i == 0:
        continue

    island = int(row[0])
    source = int(row[1])
    ra = float(row[6])
    dec = float(row[8])
    a = float(row[14]) / 3600
    b = float(row[16]) / 3600
    pa = float(row[18])

    if island not in islands:
        islands[island] = {'sources': [], 'stats': {}}

    islands[island]['sources'].append({
        'ra': ra,
        'dec': dec,
        'a': a,
        'b': b,
        'pa': pa,
    })


for key, island in islands.items():
    folder = '{}/{}'.format(FOLDER, key)
    p = pathlib.Path(folder)
    if not p.exists():
        continue

    for freq in p.iterdir():
        f = p.joinpath(freq)
        if not f.is_dir():
            continue

        for fs in f.iterdir():
            if fs.name not in ['.DS_Store', 'selected.fits']:
                fitsfile = pathlib.Path(fs)

        # print(fitsfile)

        hdulist = fits.open(fitsfile)
        wcshelper = WCSHelper(hdulist[0].header)
        image = np.squeeze(hdulist[0].data)
        mask = np.zeros(image.shape)

        for source in island['sources']:
            x, y, sx, sy, pa = wcshelper.sky2pix_ellipse(
                (source['ra'], source['dec']),
                source['a'],
                source['b'],
                source['pa']
            )
            draw_ellipse(mask, x, y, sx, sy, pa)

        # Rudimentary source detection
        bg = image.copy().ravel()
        bg = bg[~np.isnan(bg)]
        total = len(bg)

        for _ in range(3):
            mean = np.mean(bg)
            clip = 3 * np.std(bg)

            bg = bg[bg < mean + clip]
            bg = bg[bg > mean - clip]

        selected = len(bg)
        # print("Calculating bg from {}/{} pixels".format(selected, total))

        mean = np.mean(bg)
        rms = np.std(bg)

        rms_map = (image - mean) / rms
        source_median = np.median(rms_map[mask == 1])
        if source_median > 1.5:
            print("Source exists for island {} at freq {}".format(key, freq))
            print(mean, rms)
            source = image[mask == 1]
            print(np.median(source), source_median)
            island['stats'][freq] = np.sum(source)

        image[mask == 1] = 1
        hdulist[0].writeto(f.joinpath('selected.fits'), clobber=True)

        hdulist.close()









