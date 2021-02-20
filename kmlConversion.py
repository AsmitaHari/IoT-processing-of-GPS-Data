import csv
import simplekml


def decdeg2dms(dd, direction):
    if direction in ["E", "W"]:
        degree = float(dd[:3])
        minutes = float(dd[3:])
    else:
        degree = float(dd[:2])
        minutes = float(dd[2:])
    if direction in ["S", "W"]:
        return -(degree + round(minutes/60, 6))
    else:
        return degree + round(minutes/60, 6)


def makeKML(file):
    inputfile = csv.reader(open(file, 'r'))
    kml = simplekml.Kml()
    ls = kml.newlinestring(name="Journey path", tessellate=1, extrude=1,
                           altitudemode='absolute')
    ls.style.linestyle.color = simplekml.Color.hexa("Af00ffff")
    ls.style.linestyle.width = 6
    ls.style.polystyle.color = simplekml.Color.hexa("7f00ff00")

    for row in inputfile:
        if len(row) > 1 and row[0] == '$GPRMC':
            # try:
            ls.coords.addcoordinates([(decdeg2dms(row[5], row[6]),
                                       decdeg2dms(row[3], row[4]),
                                       float(row[7]))])
            # except BaseException:
            #     kml.newpoint(name=row[0], coords=[(row[3], row[1])])
    kml.save('carplaces_new.kml')
    print("Done with KML conversion")


#what data should we be including
"""def edgeCases(previousPoint,point):
    if point[]
def mapKML():
def findStops():
def findturns():
"""


def main():
    makeKML("files/ZJ42_EC0_to_RIT.TXT")


if __name__ == '__main__':
    main()