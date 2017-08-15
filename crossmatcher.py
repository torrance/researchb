#! /usr/bin/env python3

import csv
import argparse
from pathlib import Path

import numpy as np

from helpers import gcd

# Look at setting this to half a beam size at the
# respective bandwidth
DISTANCE = 0.025  # degrees


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, nargs='+')
    args = parser.parse_args()

    observations = []

    # First collect *all* observations
    bands = set()
    for file in args.file:
        path = Path(file)
        band = path.parts[-2]
        bands.add(band)

        with path.open() as f:
            print("Collecting sources from {}".format(str(path)))
            reader = csv.reader(f)
            for i, row in enumerate(reader):
                # Skip header
                if i == 0:
                    continue

                ra = float(row[6])
                ra_err = float(row[7])
                dec = float(row[8])
                dec_err = float(row[9])
                peak = float(row[10])
                total = float(row[12])
                a = float(row[14]) / 3600
                b = float(row[16]) / 3600
                wpa = float(row[18])

                obs = Observation(ra, dec, ra_err, dec_err, a, b, wpa, peak, total, band)
                observations.append(obs)

    print("Total observations: {}".format(len(observations)))

    # Now we attempt to match sources
    # Todo: no attempt made here to match sources that are transitively connected
    sources = []
    for obs in observations:
        for source in sources:
            if source.distance(obs) < DISTANCE:
                source.add(obs)
                break
        else:
            source = Source(obs)
            sources.append(source)

    print("Total sources (after merging): {}".format(len(sources)))

    # Output CSV
    with open('matched.csv', 'w') as f:
        writer = csv.writer(f)
        bands = list(bands)
        header = ['ra', 'dec', 'a', 'b', 'wpa'] + bands
        writer.writerow(header)

        for source in sources:
            row = [
                source.ra(),
                source.dec(),
                source.a(),
                source.b(),
                source.wpa()
            ]
            for band in bands:
                totals = source.totals()
                row.append(totals.get(band, ''))

            writer.writerow(row)

    # Output annotation
    with open('matched.ann', 'w') as f:
        colors = ['RED', 'BLUE', 'YELLOW', 'ORANGE', 'GREEN']
        for i, source in enumerate(sources):
            color = colors[i % 5]

            for obs in source.observations:
                print("COLOR {}".format(color), file=f)
                print("TEXT W {} {} {}".format(
                    obs.ra,
                    obs.dec,
                    obs.band
                ), file=f)
                print("ELLIPSE W {} {} {} {} {}".format(
                    obs.ra,
                    obs.dec,
                    obs.a,
                    obs.b,
                    -obs.wpa
                ), file=f)


class Source:
    def __init__(self, observation):
        self.observations = [observation]

    def add(self, observation):
        # band = observation.band
        # for i, obs in enumerate(self.observations):
        #     if obs.band == band:
        #         # Delete existing one if it is worse
        #         if observation.ra_err * observation.dec_err < obs.ra_err * obs.dec_err:
        #             del self.observations[i]
        #         else:
        #             # If the new one is worse, just return
        #             return

        self.observations.append(observation)

    def distance(self, observation):
        dist = np.inf
        for obs in self.observations:
            d = gcd(obs.ra, obs.dec, observation.ra, observation.dec)
            if d < dist:
                dist = d

        return dist

    def best_obs(self):
        # Go for highest, widest bandwidth
        best = None
        best_band = None
        for obs in self.observations:
            band = obs.band.split('-')
            band = [int(band[0]), int(band[1])]

            if best:
                # Prefer widest band always
                if (band[1] - band[0]) > (best_band[1] - best_band[0]):
                    print(band, "better than ", best_band)
                    best = obs
                    best_band = band
                # For bands of similar width, prefer highest frequency
                elif (best_band[1] - best_band[0]) - (band[1] - band[0]) < 5 and band[0] > best_band[0]:
                    print(band, "better than ", best_band)
                    best = obs
                    best_band = band
                else:
                    continue
            else:
                # Populate initial
                best = obs
                best_band = band

        return best

    def ra(self):
        ras = []
        for obs in self.observations:
            ras.append(obs.ra)

        return np.mean(ras)

    def dec(self):
        decs = []
        for obs in self.observations:
            decs.append(obs.dec)

        return np.mean(decs)

    def a(self):
        return self.best_obs().a

    def b(self):
        return self.best_obs().b

    def wpa(self):
        return self.best_obs().wpa

    def totals(self):
        totals = {}
        for obs in self.observations:
            totals[obs.band] = obs.total

        return totals


class Observation:
    def __init__(self, ra, dec, ra_err, dec_err, a, b, wpa, peak, total, band):
        self.ra = ra
        self.dec = dec
        self.ra_err = ra_err
        self.dec_err = dec_err
        self.a = a
        self.b = b
        self.wpa = wpa
        self.peak = peak
        self.total = total
        self.band = band


if __name__ == '__main__':
    main()
