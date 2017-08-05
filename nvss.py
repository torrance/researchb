#! /usr/bin/env python3

f = open('../Catalog/nvss.txt')
annf = open('../Output/nvss.ann', 'w')
print("COLOR PINK", file=annf)

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

    peak = float(parts[7])
    a = float(parts[8]) / 3600
    b = float(parts[9]) / 3600
    wpa = float(parts[10])

    print("ELLIPSE W {ra} {dec} {maj} {min} {pa}".format(
        ra=ra,
        dec=dec,
        min=b,
        maj=a,
        pa=wpa
    ), file=annf)
