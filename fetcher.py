#! /usr/bin/env python3

import argparse
from pathlib import Path

import requests


def main():
    # parser = argparse.ArgumentParser(description='Retrieves tiles from GELAM')
    # parser.add_argument('ra', type=float)
    # parser.add_argument('dec', type=float)
    # args = parser.parse_args()
    ras = [22.5, 17.5, 12.5, 7.5, 2.5, -357.5, -352.5, -347.5, -342.5, -337.5]
    decs = [-7, -12, -17, -22, -27, -32, -37, -42, -47]

    for ra in ras:
        for dec in decs:
            dl(ra, dec)


def dl(ra, dec):
    print("Processing: RA: {:.1f}, DEC: {:.1f}".format(ra, dec))
    r = requests.post("http://mwa-web.icrar.org/gleam_postage/q/form", data={
        '_charset_': 'UTF-8',
        '__nevow_form__': 'genForm',
        'pos': "{:.2f}, {:.2f};".format(ra, dec),
        'size': '5',
        'proj_opt': 'ZEA',
        '_DBOPTIONS_ORDER': 'freq',
        'MAXREC': '100',
        '_FORMAT': 'JSON',
        '_VERB': 'H',
        'submit': 'Go',
    })

    try:
        images = r.json()['data']
    except Exception as e:
        print("Failed to parse json {}".format(str(e)))
        print(r.content)
        exit()

    for freq, url in images:
        folder = Path('.').joinpath('{:.1f}-{:.1f}'.format(ra, dec), freq)
        if folder.exists():
            continue

        print("Downloading to {}".format(folder))
        r = requests.get(url)

        # Get filename
        cds = r.headers['Content-disposition'].split(';')
        for cd in cds:
            if cd[0:9] == 'filename=':
                filename = cd[9:]
                break

        if filename is None:
            filename = url

        folder.mkdir(exist_ok=True, parents=True)

        f = folder.joinpath(filename).open(mode='wb')
        f.write(r.content)
        f.close()


if __name__ == '__main__':
    main()