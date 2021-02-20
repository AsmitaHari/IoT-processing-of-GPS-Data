"""
Project- Big Data Analytics (CSCI 720)
This program process the gps data and identifies the stop,left and right turn
and uses agglomeration to combine all the left,stop and the right turns
Author - Pawan Shetty
Author - Luke Batchelder
Author- Asmita Hari
"""

import simplekml
import math
import os
import glob
from sklearn.cluster import AgglomerativeClustering
import pandas as pd, numpy as np, matplotlib.pyplot as plt, time
from sklearn.cluster import DBSCAN
from shapely.geometry import MultiPoint
from sklearn.metrics.pairwise import haversine_distances
from math import radians

percentVotes=.0

"""
Hyperparameters used for testing different stop time, angle of turn
max and min change of lat and long
"""
# Keeping all global variables here and in capitals
ORIENTATIONS = [0, 1, 2, 3]
ORIENTATION_NOTATIONS = ["west", "north", "east", "south"]
MAX_CHANGE = .5
MIN_CHANGE = 0.000005
TURN_SIZE = 90
PREV_POINTS_SIZE = 3
SKIP_QUANTITY = 1
SAME_SPEED_BUFFER = 5
STOP_VELOCITY = 0.005
MIN_STOP_TIME = 60
MAX_STOP_TIME = 120
ANGLE_BUFFER = 40

"""
Class point is used to store point read from the gps text file
"""
class Point:
    """
     Point Class, can be used to store a point read from the txt files
    """
    __slots__ = "longitude", "latitude", "speed", "angle", "time"

    def __init__(self, longitude, latitude, speed, angle, time):
        """
        constructor for the Point class can store details for each point
        :param longitude: the longitude at that point
        :param latitude: the latitude at that point
        :param speed: the speed at that point
        :param angle: the angle to the magnetic north at that point
        """
        self.longitude = longitude
        self.latitude = latitude
        self.speed = speed
        self.angle = angle
        self.time = time

"""
This class is used to store previous points for comparison
"""
class PreviousPoints:
    __slots__ = "points", "sameSpeedTacker", "skipCondition", "currentSkip"

    def __init__(self):
        self.points = [Point(0.0, 0.0, 0.0, 0.0, 0.0)] * PREV_POINTS_SIZE
        self.skipCondition = False
        self.currentSkip = SKIP_QUANTITY
    """
     This method is used to check if the speed for the given point with all other points 
    """
    def checkSameVelocity(self, newPoint):
        for x in self.points:
            if (x.speed - SAME_SPEED_BUFFER > newPoint.speed) \
                    or (x.speed + SAME_SPEED_BUFFER < newPoint.speed)\
                    or x.angle < ANGLE_BUFFER:
                return False
        return True
    """
    This method add a new point to a list of point and checks if the speeds are same
    :param newPoint: the new point to be added
    """
    def addPoint(self, newPoint):
        self.points.pop(0)
        self.points.append(newPoint)
        if self.checkSameVelocity(newPoint):
            if self.currentSkip > 0:
                self.currentSkip -= 1
                return False
            else:
                self.currentSkip = SKIP_QUANTITY
                return False
        return True

"""
This method takes the raw lat and longitude values from the 
text file calculates the degree and minutes and get the actual lat and long values
:param dd: contains the lat and long values 
:param direction: The direction of the point
:param prevPoint: contains the previous point
:return the actual lat and long values
"""
def decdeg2dms(dd, direction,prevPoint):
    try:
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
    except ValueError:
        return prevPoint.angle
"""

Function which reads the gps data for a one file mentioned and 
does the preprocessing step and gives a coords list

:param fileName: Name of the text file to be read
:return the coordsList which is the list of the coordinates

"""

def readGPSDataForOneFile(fileName):
    file = open(fileName)
    speedTracker = PreviousPoints()
    coordsList = []
    speed, angle, old_longitude, old_latitude = 0.0, 0.0, 0.0, 0.0
    lastPoint = Point(0.0, 0.0, 0.0, 0.0, 0.0)
    try:
        for line in file.readlines():
            try:
                line_data = line.split(",")
                if line_data[0] != "$GPGGA" and line_data[0] != "$GPRMC":
                    continue
                if line_data[0] == "$GPGGA":
                    longitude, latitude, time, _ = get_GGPA_data(line_data,
                                                                        speedTracker)
                else:
                    longitude, latitude, speed, angle, time, _ = \
                        get_GPRMC_data(line_data, speedTracker)
                current_point = Point(longitude, latitude, speed, angle, time)
                if lastPoint is None:
                    lastPoint = current_point

                if checkZero(lastPoint, current_point):
                    coordsList.append(current_point)
                lastPoint = current_point
            except ValueError:
                continue
    except UnicodeEncodeError and UnicodeDecodeError:
        pass
    return coordsList

"""

Function which read the gps data and finds the points which are left, right turns and stop
:param fileName: Name of the text file
:return coordsList: List of all coordinates
:return leftTurnPointsList: List of all left turns
:return rightTurnPointsList: List of all right turns 
:return stopPointsList: List of stopPoints
:return coordsRefList: List of coordinate refernce points

"""
def readGPSData(fileName):
    file = open(fileName)
    speedTracker = PreviousPoints()
    coordsRefList = {}
    coordsList, stopPointsList, leftTurnPointsList, rightTurnPointsList = [], \
                                                                          [], \
                                                                          [], \
                                                                          []
    speed, angle, old_longitude, old_latitude = 0.0, 0.0, 0.0, 0.0
    stopped = False
    lastPoint = Point(0.0, 0.0, 0.0, 0.0, 0.0)
    try:
        for line in file.readlines():
            pass
            try:
                line_data = line.split(",")
                if line_data[0] != "$GPGGA" and line_data[0] != "$GPRMC":
                    continue
                # msg = pynmea2.parse(line, check=True)
                # if msg.sentence_type == "GGA" or msg.sentence_type == "RMC":
                if line_data[0] == "$GPGGA":
                    longitude, latitude, time, turnData = get_GGPA_data(line_data,
                                                                        speedTracker)
                else:
                    longitude, latitude, speed, angle, time, turnData = \
                        get_GPRMC_data(line_data, speedTracker)
                # Creating the current point
                current_point = Point(longitude, latitude, speed, angle, time)
                if lastPoint is None:
                    lastPoint = current_point
                """
                checks two lat and long points for biggest jump and if the points are the same 
                """
                if checkZero(lastPoint, current_point):
                    """
                    Checks if the point is a stop point or not
                    """
                    if not stopped:
                        stopped = checkIfStopped(current_point, stopped, stopPointsList)
                        if stopped:
                            stopPointsList.append(current_point)
                    else:
                        stopped = checkIfStopped(current_point, stopped, stopPointsList)
                        if stopped:
                            continue

                    # checks if the turn is right or left and adds to list accordingly
                    if turnData == 'right':
                        rightTurnPointsList.append(current_point)
                        coordsRefList[(longitude, latitude)] = turnData
                    elif turnData == "left":
                        leftTurnPointsList.append(current_point)
                        coordsRefList[(longitude, latitude)] = turnData
                    """
                    if at a certain velocity dont add to to the list
                    """
                    if not speedTracker.addPoint(current_point):
                        coordsList.append(current_point)

                    speedTracker.addPoint(current_point)
                lastPoint = current_point
            except ValueError:
                continue
    except UnicodeEncodeError and UnicodeDecodeError:
        pass
    if len(stopPointsList) > 1:
        stopPointsList.pop()
    return coordsList, leftTurnPointsList, rightTurnPointsList, \
           stopPointsList, coordsRefList

"""
Function which checks if the current point's speed is same as the previous speed and checks if the 
time of stop points is less that minimum stop time and greater than the max stop time 
:param current_point: Current point being checked
:pram stopped: boolean value 
:param stopPointsList: List of stop points
:return boolean

"""
def checkIfStopped(current_point, stopped, stopPointsList):
    if current_point.speed < STOP_VELOCITY:
        if not stopped:
            return True
        else:
            return False
    else:
        if len(stopPointsList) > 1:
            prev_stop_point = stopPointsList[-1]
            if current_point.time - prev_stop_point.time < MIN_STOP_TIME \
            or current_point.time - prev_stop_point.time > MAX_STOP_TIME:
                if stopped:
                    return True
                else:
                    return False
        else:
            return True

"""
Function which reads the line from a GPRMC calcautes the degrees and mintues and finds the speed and angle
:param line_data: line from the text file
:param speedTracker : of class points used to tracking speed
:return longitude, latitude, speed, angle, time, turnData
"""
def get_GPRMC_data(line_data, speedTracker):
    longitude = round(
        decdeg2dms(line_data[5], line_data[6], speedTracker.points[-1]), 6)
    latitude = round(
        decdeg2dms(line_data[3], line_data[4], speedTracker.points[-1]), 6)
    speed = float(line_data[7])
    angle = float(line_data[8])
    a_point = speedTracker.points[-1]
    turnData = checkTurn(a_point.angle, angle, 12)
    time = float(line_data[1])
    return longitude, latitude, speed, angle, time, turnData

"""
Function which reads the line from a GGA calculates the degrees and mintues for longitude, latitude
:param line_data: line from the text file
:param speedTracker : of class points used to tracking speed 
:return longitude, latitude, time, turnData
"""
def get_GGPA_data(line_data, speedTracker):
    longitude = round(
        decdeg2dms(line_data[4], line_data[5], speedTracker.points[-1]), 6)
    latitude = round(
        decdeg2dms(line_data[2], line_data[3], speedTracker.points[-1]), 6)
    turnData = "no"
    time = float(line_data[1])
    return longitude, latitude, time, turnData

""""
Function which checks given two angles if the turn made is left turn or right turn
:param angle1: angle of point 1
:param angle2: angle of point 2
:param min_angle: hyperParamenter to check the angle difference
:return right if the turn is right ,left if the turn is left else no

"""
def checkTurn(angle1, angle2, min_angle):
    if angle2 > angle1:
        if abs(angle2 - angle1) < 180:
            if abs(angle2 - angle1) > min_angle:
                turnData = "right"
            else:
                turnData = "no"
        else:
            if abs(angle2 - angle1) > min_angle:
                turnData = "left"
            else:
                turnData = "no"
    else:
        if abs(angle1-angle2) < 180:
            if abs(angle1 - angle2) > min_angle:
                turnData = "left"
            else:
                turnData = "no"
        else:
            if abs(angle1 - angle2) > min_angle:
                turnData = "right"
            else:
                turnData = "no"
    return turnData
""""
This methods adds the lat and long values from all the points to a final list
:param: pinLists: List of all the pins
:return finalList: List of Lat and long values
"""

def buildListFromDict(pinLists, count, percentVotes):
    finalList = []
    for x in pinLists.keys():
        # if pinLists[x]/count>=percentVotes:
        finalList.append([x[0], x[1]])
    return finalList
"""
Function which creates a final Kml file by reading in all the coordslist,
stops,lefts and rights and assigns  color to each points and a pin in the kml file

:param coordsList: List of all the coordinates
:param stops: List of all the stop points
:param leftTurn: Lsit of all the left turns
:param rightTurns: List of all the right turns 
:param name: Name of the kml file

"""

def getKML(coordsList, stops, leftTurns, rightTurns, name):

    kml = simplekml.Kml()
    stops = agglomerateLists(stops) ## calls the agglomeration on the stop points
    leftTurns = agglomerateLists(leftTurns) # calls the agglomeration on the right points
    rightTurns = agglomerateLists(rightTurns) # calls the agglomeration on the left points
    total_len = len(stops) + len(rightTurns) + len(leftTurns)
    print("Total no of points:", total_len)
    for x in stops: # for all the stop points names the pin as stop and assigns  red color
        pnt = kml.newpoint(name="stop", coords=[x])
        pnt.style.iconstyle.color = '64780078'
        pnt.style.labelstyle.color = '64780078'  # Red
    for x in leftTurns: # for all the left points names the pin as left and assigns  yellow color
        pnt = kml.newpoint(name="LEFT", coords=[x])
        pnt.style.iconstyle.color = '6414F0FF'
        pnt.style.labelstyle.color = '6414F0FF'
    for x in rightTurns: # for all the right points names the pin as left and assigns  yellow color
        pnt = kml.newpoint(name="RIGHT", coords=[x])
        pnt.style.iconstyle.color = '64F0FF14'
        # why cyan for right turn and paths?
        pnt.style.labelstyle.color = '64F0FF14'
    ls = kml.newlinestring(name="Description", tessellate=1,
                           extrude=1, altitudemode='absolute')
    ls.style.linestyle.color = simplekml.Color.hexa("14FFF0FF")
    ls.style.linestyle.width = 6
    ls.style.polystyle.color = simplekml.Color.hexa("00F0145A")
     ## adds in all the coordinates
    for k in coordsList:
        ls.coords.addcoordinates(coords=[get_formatted_coords(x) for x
                                         in k])
    kml.save("kml/"+name+".kml")
    print("Done with KML for " + name)


def dbscan(a_list):
    kms_per_radian = 6371.0088
    df = pd.DataFrame(a_list)
    coords = df[[0, 1]].to_numpy()
    epsilon = 1.5 / kms_per_radian

    db = DBSCAN(eps=epsilon, min_samples=10, algorithm='ball_tree',
                metric='haversine').fit(np.radians(coords))
    labels_list = db.labels_

    # get the number of clusters
    # num_clusters = len(set(cluster_labels))
    #
    # # all done, print the outcome
    # message = 'Clustered {:,} points down to {:,} clusters, for {:.1f}% compression in {:,.2f} seconds'
    # print(message.format(len(df), num_clusters, 100 * (1 - float(num_clusters) / len(df)), time.time() - start_time))
    # print('Silhouette coefficient: {:0.03f}'.format(metrics.silhouette_score(coords, cluster_labels)))
    #
    # clusters = pd.Series([coords[cluster_labels == n] for n in range(num_clusters)])
    # return centerCal(clusters, df)
    return get_centroids(a_list, labels_list)


def centerCal(clusters, df):
    centermost_points = clusters.map(get_centermost_point)
    lats, lons = zip(*centermost_points)

    # from these lats/lons create a new df of one representative point for each cluster
    rep_points = pd.DataFrame({'lat': lats, 'lon': lons})
    return rep_points

def get_centermost_point(cluster):
        centroid = (MultiPoint(cluster).centroid.x, MultiPoint(cluster).centroid.y)
"""
Function which uses the sklearn agglomeration and returns the a new list of points which are merged together using
the centriod
:param a_list: list of latitude and longitude values
:return a new list of points 
"""
def agglomerateLists(a_list):
    an_array = np.array(a_list)
    clustering = AgglomerativeClustering(
        n_clusters=None, distance_threshold=0.0009).fit(
        an_array)
    labels_list = clustering.labels_
    return get_centroids(a_list, labels_list)

""""
Function which finds the centroid between the lat and long values
:param a_list: list of latitude and longitude values
:param labels_list : list of labels from the agglomeration clustering
"""
def get_centroids(a_list, labels_list):
    a_dict = {}
    centroid_list = []
    for index, a_cluster in np.ndenumerate(labels_list):
        index = index[0]
        if a_cluster not in a_dict:
            a_dict[a_cluster] = [a_list[index]]
        else:
            temp_list = a_dict[a_cluster]
            temp_list.append(a_list[index])
            a_dict[a_cluster] = temp_list
    for a_key in a_dict:
        points_list = a_dict[a_key]
        temp_list = sorted(points_list, key=lambda x: (math.sqrt(x[0]**2 + x[
            1]**2)))
        if len(temp_list) > 2:
            centroid_list.append(temp_list[len(temp_list)//2])
    return centroid_list

"""
function which reads the point from the point class and returns a tuple 
of longitude,latitude and speed
:param a_point: A point from the point class
"""
def get_formatted_coords(a_point):
    return [a_point.longitude, a_point.latitude, a_point.speed]


"""
Fucntion which pre process the given lat and long values and check if the two are same
and if there is huge jump in the values

:param point1: point class which contains lat,long speed etc
:param point2: point class which contains lat,long speed etc

:return boolean
"""
def checkZero(point1, point2):
    abs_longitude = abs(point1.longitude - point2.longitude)
    abs_latitude = abs(point1.latitude - point2.latitude)
    if (abs_latitude + abs_longitude) ** 2 <= MAX_CHANGE \
           and \
           ((abs_latitude + abs_longitude) > MIN_CHANGE):
        return True
    else:
        return False
"""

Function which checks if the given tuple point exits then increments the count of the element in the dictionary
:param dict: dictionary of tuples
:param list: list of points
 
"""

def addToDict(dict, list):
    for x in list:
        hashableElement = tuple(get_formatted_coords(x))
        if hashableElement in dict:
            dict[hashableElement] += 1
        else:
            dict[hashableElement] = 0

"""
Function which creates a kml foler and asks the user to input the filename of the .txt file and reads in the  
gps data and gives a kml file called Fina;PathFile
"""
def getKMLPath():
    if not os.path.isdir('kml'):
        os.mkdir('kml')
    fullCoordsList, fullcoordsRefenceDict = [], []
    file_name = input("Enter the filename:")
    for filePath in glob.glob(file_name):
        coordsList = readGPSDataForOneFile(filePath)
        fullCoordsList.append(coordsList)
    getKML(fullCoordsList, [], [], [], "FinalPathFile")

"""

This method :
creates a folder called kml for all the kml files
Reads all the .txt files from the folder
Finds all the points in the txt file which are left ,right turns and stop points
combines all the points and calls the getkml function to generate a final kml file
"""

def getStopsandTurns():
    count = 0
    if not os.path.isdir('kml'):
        os.mkdir('kml')
    fullCoordsList, fullLeftTurnPointsList, fullrightTurnPointsList, \
        fullstopPointsList, fullcoordsRefenceDict = [], {}, {}, {}, []
    directory_name = input("Enter the path to the directory were all files "
                           "are stored:")
    directory_name = directory_name.replace("\"", "") + "/*.txt"
    for filePath in glob.glob(directory_name):
        count += 1
        coordsList, leftTurnPointsList, rightTurnPointsList, stopPointsList, coordsRefenceList = \
            readGPSData(filePath)
        addToDict(fullLeftTurnPointsList, leftTurnPointsList)
        addToDict(fullrightTurnPointsList, rightTurnPointsList)
        addToDict(fullstopPointsList, stopPointsList)
    fullLeftTurnPointsList = buildListFromDict(fullLeftTurnPointsList, count,
                                               percentVotes)
    fullrightTurnPointsList = buildListFromDict(fullrightTurnPointsList, count,
                                                percentVotes)
    fullstopPointsList = buildListFromDict(fullstopPointsList, count, percentVotes)
    getKML([], fullLeftTurnPointsList, fullrightTurnPointsList,
           fullstopPointsList, "FinalStopsAndTurnsFile")
"""

Main method asks the user for input 
choice a if the user wants to Convert one GPS file to one KML map showing the path traveled
choice b if the user wants to Convert a set of GPS files to one KML map of stops and turns

"""

def main():
    choice = "None"
    while choice == "None":
        print("What would you like to do (Type a or b based on your choice)?")
        print("a) Convert one GPS file to one KML map showing the path traveled")
        print("b) Convert a set of GPS files to one KML map of stops and turns")
        print("")
        choice = input("Enter your choice")
        if choice == "a":
            getKMLPath()
        elif choice == "b":
            getStopsandTurns()
        else:
            print("Wrong choice, Try Again!")
            choice = "None"

    print("************ Done generating KML files ***************")


if __name__ == '__main__':
    main()