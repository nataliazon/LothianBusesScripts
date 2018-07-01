from collections import defaultdict
import json
import utm
import re
import math
from string import ascii_uppercase
import itertools

import re


##########
# README
# Look below for the main() function to include the path to the files with the services and stops.
# These are Json formatted txt files, obtained from Transport For Edinburgh Data API
#########

def convertToEastingNorthing(latitude, longitude):
    converted = (utm.from_latlon(latitude, longitude))
    easting = converted[0]
    northing = converted[1]
    return([easting, northing])


def loadJsonFile(pathToFile):
    fullFile = open(pathToFile).read()

    beg = fullFile.find("<")
    end = fullFile.rfind(">")

    correctJsonFile = fullFile[end + 1:]

    loadedFileData = json.loads(correctJsonFile)
    return loadedFileData


def getRoutePointsByServiceByDestination(pathToServiceFile):
    data = loadJsonFile(pathToServiceFile)
    listOfServiceData = data['services']
    serviceNames = set()
    serviceDataByNames = dict()
    serviceDestinationsByServiceNumber = dict()
    for serviceItem in listOfServiceData:
        serviceDestinations = set()
        serviceNames.add(serviceItem['name'])
        serviceDataByNames[serviceItem['name']] = serviceItem

        for routeData in serviceItem['routes']:
            serviceDestinations.add(routeData['destination'])
        serviceDestinationsByServiceNumber[serviceItem['name']] = serviceDestinations



    routePointsByServiceByDestination = dict()

    for serviceName in serviceDataByNames.keys():
        routesData = serviceDataByNames[serviceName]['routes']
        routePointsByServiceByDestination[serviceName] = dict()
        for destination in serviceDestinationsByServiceNumber[serviceName]:

            for routeData in routesData:
                if routeData['destination'] == destination:
                    routePoints = routeData['points']
                    routePointsByServiceByDestination[serviceName][destination] = routePoints


    return routePointsByServiceByDestination


def getConvertedRoutePointsByServiceByDestination(routesData):
    convertedByServiceByDestination = dict()
    for service in routesData.keys():
        convertedByServiceByDestination[service]=dict()
        for destination in routesData[service].keys():

            routePoints = routesData[service][destination]
            convertedRoutePoints= list()
            for pointData in routePoints:
                stopId = pointData['stop_id']
                latitude = pointData['latitude']
                longitude = pointData['longitude']
                convertedLocation = convertToEastingNorthing(latitude, longitude)
                easting = convertedLocation[0]
                northing = convertedLocation[1]
                convertedRoutePoints.append({'stop_id':stopId, 'easting':easting, 'northing':northing})
            convertedByServiceByDestination[service][destination] = convertedRoutePoints
    return convertedByServiceByDestination


def mapIdsToBusStops(routesData, stopsFilePath):
    loadedStops = loadJsonFile(stopsFilePath)
    nameOfStopById = dict()
    for item in loadedStops['stops']:
        nameOfStopById[item['stop_id']] = item['name']

    dataWithNames = dict()
    for service in routesData.keys():
        dataWithNames[service]=dict()
        for destination in routesData[service].keys():
            for item in routesData[service][destination]:
                st_id = item['stop_id']
                if st_id is not None and st_id in nameOfStopById.keys():

                    name = nameOfStopById[st_id]
                    item['name'] = name
                elif st_id is not None:
                    print("!! Service " + str(service) + " to " + str(destination) + ": stop ID " + str(st_id) + str(" not in stops.txt data. Skipping ..."))

                else:
                    item['name'] = "None"
    return routesData


def generateCarmaXYtabsAndMappingFile(data, service, destination, startID, pathToMappingFileOutput):

    if service not in data.keys():
        raise Exception("No such service " + service)
    if destination not in data[service].keys():
        raise Exception("No such destination " + destination + " for service " + service)

    requestedData = data[service][destination]

    indexValue = startID

    eastings = dict()
    northings = dict()

    stopNames = dict()
    stopIds = dict()

    for item in requestedData:
        if item['stop_id'] == 36236382:
            continue
        eastings[indexValue] = item['easting']
        northings[indexValue] = item['northing']
        stopNames[indexValue] = item['name']
        stopIds[indexValue] = item['stop_id']

        indexValue +=1

    outf = open(pathToMappingFileOutput, "w")

    stringToWrite = "Service " + service + " to " + destination + "\n"
    stringToWrite += "carma_id;api_id;name;easting;northing\n\n"
    for idv in stopNames.keys():
        stringToWrite += (str(idv) + ";" + str(stopIds[idv]) + ";" + str(stopNames[idv]) + ";" + str(eastings[idv]) + ";" + str(northings[idv])) + "\n"
        print("if( location_current == " + str(idv) + ") {")
        if(stopIds[idv] != None):
            print("    return true;}")
        else:
            print("    return false;}")

    print(stringToWrite)
    outf.write(stringToWrite)
    outf.close()

    carmaXstr = "const x = [: "
    carmaYstr = "const y = [: "

    for idv in stopNames.keys():
        carmaXstr += "{:.4f}".format(eastings[idv]) + ", "
        carmaYstr += "{:.4f}".format(northings[idv]) + ", "


    carmaXstr = carmaXstr[:-2] + " :];\n"
    carmaYstr = carmaYstr[:-2] + " :];\n"

    prefix = "// BEGIN: Auto-generated location coordinates for service " + service + " to " + destination +".\n "
    suffix = "// END: Auto-generated location coordinates. "

    return prefix + carmaXstr + carmaYstr + suffix


def main():
    print("Running locations generator script...")
    pathToServicesFile = "/Users/natalia/Work/LothianBusesData/services.txt"
    pathToStopsFile = "/Users/natalia/Work/LothianBusesData/stops.txt"


    routesData = getRoutePointsByServiceByDestination(pathToServicesFile)
    convertedRoutesData = getConvertedRoutePointsByServiceByDestination(routesData)
    dataWithNames = mapIdsToBusStops(convertedRoutesData, pathToStopsFile)

    service = "5"
    destination = "Hunter's Tryst"
    startID = 0
    pathToMappingFileOutput = "/Users/natalia/Work/mapping.txt"
    carmaXYtext = generateCarmaXYtabsAndMappingFile(dataWithNames, service, destination, startID, pathToMappingFileOutput)

    print(carmaXYtext)








if __name__ == "__main__":
    main()