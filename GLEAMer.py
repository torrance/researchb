#! /usr/bin/env python3

import csv
import pathlib

from astropy.io import fits
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

from helpers import WCSHelper, draw_ellipse


CATALOG = '/Users/torrance/Google Drive/Research B/Output/blanked_comp.csv'
FOLDER = '/Users/torrance/Google Drive/Research B/Stamps'
OUTPUT = '/Users/torrance/Google Drive/Research B/Output/'

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

freqs = set()

for key, island in islands.items():
    print("Processing {}".format(key), end="\r")
    folder = '{}/{}'.format(FOLDER, key)
    p = pathlib.Path(folder)
    if not p.exists():
        continue

    for freq in p.iterdir():
        f = p.joinpath(freq)
        if not f.is_dir():
            continue

        freqs.add(freq.name)

        for fs in f.iterdir():
            if fs.name not in ['.DS_Store', 'selected.fits']:
                fitsfile = fs

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

        source = image[mask == 1]
        island['stats'][freq.name] = {'flux': np.sum(source), 'median': source_median}

        #image[mask == 1] = 1
        #hdulist[0].writeto(f.joinpath('selected.fits'), clobber=True)

        hdulist.close()

print("")

for key, island in islands.items():
    wideband_x = []
    wideband_y = []
    narrowband_x = []
    narrowband_y = []

    for freq, stats in island['stats'].items():
        if stats['median'] < 1.5:
            continue

        low, high = freq.split('-')
        low = int(low)
        high = int(high)
        mid = (low + high) / 2

        # Wideband data
        if high-low > 15:
            wideband_x.append(mid)
            wideband_y.append(stats['flux'])
        # Narrowband data
        else:
            narrowband_x.append(mid)
            narrowband_y.append(stats['flux'])

    combined_x = wideband_x + narrowband_x
    combined_y = wideband_y + narrowband_y

    if len(combined_x) > 7:
        plt.figure('island-{}'.format(key), figsize=(8, 10))

        slope, intercept, r_value, p_value, stderr = linregress(np.log(combined_x), np.log(combined_y))
        island['spectral_index'] = {'slope': slope, 'stderr': stderr, 'rsquared': r_value**2, 'p_value': p_value}

        ax = plt.subplot(2, 1, 1)
        plt.scatter(wideband_x, wideband_y, color='red')
        plt.scatter(narrowband_x, narrowband_y, color='blue')
        x_low = min(combined_x)
        x_high = max(combined_x)
        plt.plot([x_low, x_high], [np.exp(intercept) * x_low**slope, np.exp(intercept) * x_high**slope])
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.xlim([x_low * 0.9, x_high * 1.1])

        plt.title('Island {}, Slope {:.2f}+\-{:.2f}, rsquared {:.2f}, pvalue: {:.2f}'.format(key, slope, stderr, r_value**2, p_value))
        plt.xlabel('Freqency (MHz)')
        plt.ylabel('Total flux (Jy)')

        # Plot residual
        wideband_residual = np.log(np.array(wideband_y)) - (np.log(np.array(wideband_x)) * slope + intercept)
        narrowband_residual = np.log(np.array(narrowband_y)) - (np.log(np.array(narrowband_x)) * slope + intercept)
        ax = plt.subplot(2, 1, 2)
        plt.scatter(np.log(wideband_x), wideband_residual, color='red')
        plt.scatter(np.log(narrowband_x), narrowband_residual, color='blue')
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.xlim([np.log(x_low * 0.9), np.log(x_high * 1.1)])
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')

        plt.title("Residual from linear regression")
        plt.xlabel('Log frequency (MHz)')
        plt.ylabel('Log residual total flux (Jy)')

        p = pathlib.Path(OUTPUT).joinpath('plots')
        p.mkdir(exist_ok=True, parents=True)
        p = p.joinpath('island-{}.pdf'.format(key))

        plt.tight_layout()
        plt.savefig(str(p))
        plt.close()
    else:
        island['spectral_index'] = {'slope': '', 'stderr': '', 'rsquared': '', 'p_value': ''}

    if len(island['sources']) > 1:
        island['type'] = 'multiple'
    elif (island['sources'][0]['a'] * island['sources'][0]['b']) / (148.3 / 3600)**2 >= 1.1:
        island['type'] = 'extended'
    else:
        island['type'] = 'single'

with pathlib.Path(OUTPUT).joinpath('gleamer.csv').open('w') as f:
    writer = csv.writer(f)

    freqs = list(freqs)  # Cast to list to 'lock in' order
    cols = ['island', 'type', 'spectral_slope', 'spectral_stderr', 'spectral_rsquared', 'spectral_pvalue']
    for freq in freqs:
        cols.append(freq)
        cols.append("{} (median RMS)".format(freq))
    writer.writerow(cols)

    for key, island in islands.items():
        row = [
            key,
            island['type'],
            island['spectral_index']['slope'],
            island['spectral_index']['stderr'],
            island['spectral_index']['rsquared'],
            island['spectral_index']['p_value'],
        ]
        for freq in freqs:
            freq = island['stats'].get(freq, {'flux': '', 'median': ''})
            row.append(freq['flux'])
            row.append(freq['median'])

        writer.writerow(row)









