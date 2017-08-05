#! /usr/bin/env python3

import requests
import csv
import pathlib

CATALOG = '/Users/torrance/Google Drive/Research B/Output/blanked_comp.csv'
FOLDER = '/Users/torrance/Google Drive/Research B/Stamps'

reader = csv.reader(open(CATALOG))

islands = []
for i, row in enumerate(reader):
    # Skip headers
    if i == 0:
        continue

    island = row[0]
    ra = row[4]
    dec = row[5]

    # Download one stamp per island only
    if island in islands:
        continue
    islands.append(island)

    r = requests.post("http://mwa-web.icrar.org/gleam_postage/q/form", data={
        '_charset_': 'UTF-8',
        '__nevow_form__': 'genForm',
        'pos': "{}, {};".format(ra, dec),
        'freq': ['072-103', '103-134', '139-170', '170-231'],
        'size': '1.5',
        'proj_opt': 'ZEA',
        '_DBOPTIONS_ORDER': 'freq',
        'MAXREC': '100',
        '_FORMAT': 'JSON',
        '_VERB': 'H',
        'submit': 'Go',
    })

    for freq, url in r.json()['data']:
        print(freq, url)
        r = requests.get(url)

        # Get filename
        cds = r.headers['Content-disposition'].split(';')
        for cd in cds:
            if cd[0:9] == 'filename=':
                filename = cd[9:]
                break

        if filename is None:
            filename = url

        folder = "{}/{}/{}".format(FOLDER, island, freq)
        print("Downloading to {}".format(folder))

        p = pathlib.Path(folder)
        p.mkdir(exist_ok=True, parents=True)

        p = pathlib.Path("{}/{}".format(folder, filename))
        f = p.open(mode='wb')
        f.write(r.content)

    if i == 4:
        break
