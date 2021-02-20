import csv
import simplekml
import numpy as np

inputfile = csv.reader(open("files/ZJ42_EC0_to_RIT.TXT", 'r'))
kml=simplekml.Kml()
points = []
for row in inputfile:
    if len(row) > 1:
        try:
            coord = (row[2], row[1]) # lon, lat order
            pnt = kml.newpoint(name=row[0], coords=[coord])
            points.append(coord)
            pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_square.png'
        except BaseException:
            kml.newpoint(name=row[0], coords=[(row[3], row[1])])

ls = kml.newlinestring(name='A LineString')
ls.coords = np.array(points)
ls.altitudemode = simplekml.AltitudeMode.relativetoground
ls.extrude = 1

kml.save("Points_and_Line.kml")