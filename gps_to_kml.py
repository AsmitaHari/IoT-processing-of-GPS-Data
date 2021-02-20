import pynmea2
import simplekml
import os
import glob


def readGPSData(fileName):
    file = open(fileName)
    coordsList = []
    speed = 0.0
    old_longitude, old_latitude = 0.0, 0.0
    for line in file.readlines():
        try:
            msg = pynmea2.parse(line, check=True)
            try:
                if msg.sentence_type == "GGA" or msg.sentence_type == "RMC":
                    longitude = round(msg.longitude, 6)
                    latitude = round(msg.latitude, 6)
                    if msg.sentence_type == "RMC":
                        speed = msg.spd_over_grnd
                    if checkZero(msg):
                        if old_longitude == longitude and old_latitude == latitude:
                            continue
                        if speed == 0.0:
                            continue
                        coordsList.append(["\n" + str(longitude), latitude, speed])
                        old_longitude = longitude
                        old_latitude = latitude
            except:
             pass
        except pynmea2.ParseError as e:
            continue

    return coordsList


def getKML(coordsList, name):

    kml = simplekml.Kml()
    ls = kml.newlinestring(name="Description", tessellate=1,
                           extrude=1, altitudemode='absolute')
    ls.style.linestyle.color = simplekml.Color.hexa("14FFF0FF")
    ls.style.linestyle.width = 6
    ls.style.polystyle.color = simplekml.Color.hexa("00F0145A")
    ls.coords.addcoordinates(coords=coordsList)
    kml.save("kml/"+name+".kml")
    print("Done with KML for " + name)


# pre processing stuff
def checkZero(msg):
    return msg.latitude != 0.0 and msg.longitude != 0.0


def main():
    if not os.path.isdir('kml'):
        os.mkdir('kml')

    print(len(glob.glob("files/*.txt")))
    for filePath in glob.glob("files/*.txt"):

        coordsList = readGPSData(filePath)
        getKML(coordsList, os.path.splitext(os.path.basename(filePath))[0])

    print("************ Done generating KML files ***************")

    # filePath = "files/2019_06_21__1422_30.txt"
    # coordsList = readGPSData(filePath)
    # getKML(coordsList, os.path.splitext(os.path.basename(filePath))[0])


if __name__ == '__main__':

    main()
