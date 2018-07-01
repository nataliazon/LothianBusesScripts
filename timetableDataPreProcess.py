__author__ = 'natalia'
from collections import defaultdict
import simplejson
import string
import json
import pickle
import numpy as np

def getAllStopNumbers(pathToStopsFile):
    fullFile = open(pathToStopsFile).read()
    stopsLoaded = json.loads(fullFile)
    allStops = stopsLoaded["stops"]
    allStopIds = list()
    for stopItem in allStops:
        allStopIds.append(stopItem['stop_id'])
    return allStopIds


def timeStrToMin(timeStr):
    if len(timeStr) != 5:
        raise Exception("Wrong string passed to timeStrToMin function")
    hoursStr  = timeStr[0:2]
    minutessStr  = timeStr[3:5]

    minutesTotal = int(hoursStr)*60 + int(minutessStr)
    if minutesTotal > 1440 or minutesTotal < 0:
        raise Exception ("Wrong value in timeStrToMin function")
    return minutesTotal


def getListWithSamplePointInsertedAndItsIndex(listOfDepartures, sampledTimeMin):
    listWithSampledPointInserted = list()
    sampledTimeIsOnList = False
    sampledTimeListId = None
    for time in sorted(listOfDepartures):
        currentDepartureTimeInMin = timeStrToMin(time)
        if (currentDepartureTimeInMin < sampledTimeMin):
            listWithSampledPointInserted.append(currentDepartureTimeInMin)
        if (currentDepartureTimeInMin == sampledTimeMin ):
            listWithSampledPointInserted.append(currentDepartureTimeInMin)
            sampledTimeListId = listWithSampledPointInserted.index(currentDepartureTimeInMin)
            sampledTimeIsOnList = True

        if (currentDepartureTimeInMin > sampledTimeMin):
            break


    if sampledTimeIsOnList is False:
        listWithSampledPointInserted.append(sampledTimeMin)
        sampledTimeListId = listWithSampledPointInserted.index(sampledTimeMin)

    for time in sorted(listOfDepartures):
        currentDepartureTimeInMin = timeStrToMin(time)
        if (currentDepartureTimeInMin > sampledTimeMin):
            listWithSampledPointInserted.append(currentDepartureTimeInMin)


    if(not sampledTimeIsOnList and not len(listWithSampledPointInserted) == len(listOfDepartures)+1):
        raise Exception("Wrong length of list")
    if(sampledTimeIsOnList and not len(listWithSampledPointInserted) == len(listOfDepartures)):
        raise Exception("Wrong length of list")
    if(sampledTimeListId == None):
        raise Exception("Wrong value of sampled point Id")

    return [listWithSampledPointInserted, sampledTimeListId, sampledTimeIsOnList]



def getFrequencyAsAverageFromPeriodOfTime(listOfDepartures,periodToAverage, sampledTimeMin):

    listWithSampledPointInfo = getListWithSamplePointInsertedAndItsIndex(listOfDepartures, sampledTimeMin)
    listWithSampledPoint = listWithSampledPointInfo[0]
    idOfSampledPoint = listWithSampledPointInfo[1]
    sampledPointWasAlreadyOnList = listWithSampledPointInfo[2]

    valuesToAverageList = list()
    if not sampledPointWasAlreadyOnList:
        for i in range(1, idOfSampledPoint -1 , 1):
            valA = float(listWithSampledPoint[idOfSampledPoint - i])
            valB = float(listWithSampledPoint[idOfSampledPoint - i -1])

            averageAB = (valA + valB)/2.0

            if averageAB > (float(listWithSampledPoint[idOfSampledPoint]) - float (periodToAverage/2.0)):
                differenceVal =  valA - valB
                valuesToAverageList.append(differenceVal)

        for i in range(idOfSampledPoint +1, len(listWithSampledPoint) -1, 1):

            valA = float(listWithSampledPoint[  i +1])
            valB = float(listWithSampledPoint[  i ])

            averageAB = (valA + valB)/2.0

            if averageAB < (float(listWithSampledPoint[idOfSampledPoint]) + float (periodToAverage/2.0)):
                differenceVal =  valA - valB
                valuesToAverageList.append(differenceVal)
    elif sampledPointWasAlreadyOnList:
        for i in range(0, idOfSampledPoint , 1):
            valA = float(listWithSampledPoint[idOfSampledPoint - i])
            valB = float(listWithSampledPoint[idOfSampledPoint - i -1])

            averageAB = (valA + valB)/2.0

            if averageAB > (float(listWithSampledPoint[idOfSampledPoint]) - float (periodToAverage/2.0)):
                differenceVal =  valA - valB
                valuesToAverageList.append(differenceVal)

        for i in range(idOfSampledPoint, len(listWithSampledPoint) -1, 1):

            valA = float(listWithSampledPoint[  i +1])
            valB = float(listWithSampledPoint[  i ])

            averageAB = (valA + valB)/2.0

            if averageAB < (float(listWithSampledPoint[idOfSampledPoint]) + float (periodToAverage/2.0)):
                differenceVal =  valA - valB
                valuesToAverageList.append(differenceVal)

    if len(valuesToAverageList) == 0:
        return 5000
    else:
        return sum(valuesToAverageList)/len(valuesToAverageList)



def getFrequencyAsWeightedAverage(listOfDepartures, sampledTimeMin):
    listWithSampledPointInfo = getListWithSamplePointInsertedAndItsIndex(listOfDepartures, sampledTimeMin)
    listWithSampledPoint = listWithSampledPointInfo[0]
    idOfSampledPoint = listWithSampledPointInfo[1]
    sampledPointWasAlreadyOnList = listWithSampledPointInfo[2]

    weightMultiplier = 1.0
    m = weightMultiplier

    w1 = lambda k,i,t_s,data : m/((abs((data[k-i] + data[k-i-1])/2.0 - t_s))**2)
    w2 = lambda k,i,t_s,data : m/((abs((data[i+1] + data[i])/2.0 - t_s))**2)

    k = idOfSampledPoint
    dat = listWithSampledPoint
    n = len(listWithSampledPoint)-1


    if k >= len(listWithSampledPoint)-1 or k <= 0:
        return 5000



    sum1 = 0.0
    for i in range(1, k-1, 1):
        sum1 += (dat[k-i] - dat[k-i-1])*w1(k, i,dat[k], dat)

    sum2 = 0.0
    for i in range(k+1,n-1, 1):
         sum2 += (dat[i+1]-dat[i])*w2(k, i,dat[k], dat)


    #print(str(k+1))
    #print(str(k-1))
    #print(str(dat))
    if not sampledPointWasAlreadyOnList:
        wK = 1.0
        itemK = abs (dat[k+1] - dat[k-1])*wK
        weightItemK = wK

    else:
        wK = 1.0
        itemK = abs(dat[k + 1] - dat[k ]) * wK + abs(dat[k] - dat[k -1 ]) * wK
        weightItemK = 2*wK


    sumW1 = 0.0
    for i in range(1, k-1, 1):
        sumW1 += w1(k, i,dat[k], dat)

    sumW2 = 0.0
    for i in range(k+1,n-1, 1):
        sumW2+= w2(k, i,dat[k], dat)

    if(sumW1 + sumW2) == 0:
        return 0.0
    result = (sum1 + sum2 + itemK)/(sumW1 + sumW2 + weightItemK)

    return result







def getFrequencyFromListOfDepartures(listOfDepartures):

    averageVal = getFrequencyAsAverageFromPeriodOfTime(listOfDepartures,3*60, 9*60)
    averageWeightedVal = getFrequencyAsWeightedAverage(listOfDepartures, 9*60)

    print(listOfDepartures)


    for i in range(0,len(listOfDepartures)-4,4):
        print(str(listOfDepartures[i]) + " & " + str(listOfDepartures[i+1]) + " & " + str(listOfDepartures[i+2]) + " & " + str(listOfDepartures[i+3]) + " \\\\ " )


    print("Average Freq at 9am with 3h period: " + str(averageVal) )
    print("Average Weighted freq at 9am : " + str(averageWeightedVal) )

    list_of_times_in_min = list()
    list_of_frequencies_avg = list()
    list_of_frequencies_weighted_avg = list()

    for hour in np.arange(0,24.5,0.25):
        minutes = hour*60.0
        aVal = getFrequencyAsAverageFromPeriodOfTime(listOfDepartures,3*60, minutes)
        aWVal = getFrequencyAsWeightedAverage(listOfDepartures, minutes)
        print( " hour " + str(hour) + " avg " + str('%.2f' % aVal) + " wAvg " + str('%.2f' % aWVal) )
        list_of_times_in_min.append(minutes)
        list_of_frequencies_avg.append(aVal)
        list_of_frequencies_weighted_avg.append(aWVal)

    return [list_of_times_in_min, list_of_frequencies_avg, list_of_frequencies_weighted_avg]

def processASingleTimetableForBusStop(number, pathPrefix, pathSuffix):
    print()
    print()
    print("Processing bus stop " + str(number))

    pathToTimetable = pathPrefix + str(number) + pathSuffix

    fullFile = open(pathToTimetable).read()

    beg = fullFile.find("<")
    end = fullFile.rfind(">")

    correctJsonFile = fullFile[end+1:]

    timatableLoaded = simplejson.loads(correctJsonFile)

    valid_From_set = set()

    departuresByService = defaultdict(list)

    for item in timatableLoaded['departures']:
        print(item['service_name'])
        departuresByService[item['service_name']].append(item)
        #print("Adding " +str(item['service_name']))
        valid_From_set.add(item['valid_from'])

    print("keys"  + departuresByService.keys())
    print("Valid from")
    print(str(valid_From_set))

    departuresByServiceByDestination = defaultdict(dict)

    setOfAllDestinations = set()

    setOfAllDayTypes = set()

    for serviceName in departuresByService.keys():
        print("#" + str(serviceName))
        listOfAllDeparturesForCurrentService = departuresByService[serviceName]
        departuresByDestination = defaultdict(list)
        for departureItem in listOfAllDeparturesForCurrentService:
            departuresByDestination[departureItem['destination']].append(departureItem)
            setOfAllDestinations.add(departureItem['destination'])
            setOfAllDayTypes.add(departureItem['day'])
        departuresByServiceByDestination[serviceName] = departuresByDestination


    departuresByServiceByDestinationByDay = defaultdict(dict)

    for serviceNo in departuresByServiceByDestination.keys():
        print("?" + str(serviceNo))
        for destinationName in departuresByServiceByDestination[serviceNo].keys():
            listOfDepartures = departuresByServiceByDestination[serviceNo][destinationName]
            departuresByDay = defaultdict(list)
            for item in listOfDepartures:
                departuresByDay[item['day']].append(item['time'])
            for item in listOfDepartures:
                departuresByDay[item['day']] = sorted(departuresByDay[item['day']])
            departuresByServiceByDestinationByDay[serviceNo][destinationName] = departuresByDay




    for service in departuresByServiceByDestinationByDay.keys():
        print(service)
        for destination in departuresByServiceByDestinationByDay[service].keys():
            for day in departuresByServiceByDestinationByDay[service][destination].keys():
                print("Working on service " + str(service) + " to " + str(destination) + " days " + str(day) + " ...")
                if(str(service) !="5" or str(destination) != "Hunter's Tryst" or str(day) != "0"):
                    print ("Not")
                    continue
                else:
                    print("YES!")
                listOfDepartures = departuresByServiceByDestinationByDay[service][destination][day]
                departuresByServiceByDestinationByDay[service][destination][day] = dict()
                departuresByServiceByDestinationByDay[service][destination][day]['departuresList'] = listOfDepartures
                result = getFrequencyFromListOfDepartures(listOfDepartures)
                departuresByServiceByDestinationByDay[service][destination][day]['frequency'] = result
                print("Times [min]: " + str(result[0]))
                print("Average delay [min]: " + str(result[1]))
                print("Weighted delay [min]: " + str(result[2]))



    # for service in departuresByServiceByDestinationByDay.keys():
    #     for destination in departuresByServiceByDestinationByDay[service].keys():
    #         for day in departuresByServiceByDestinationByDay[service][destination].keys():
    #             print(str(service) + " " +destination+ " " + str(day) + str((departuresByServiceByDestinationByDay[service][destination][day])))

    return departuresByServiceByDestinationByDay

def main(pathToStopsFile, prefixPathTimetables, suffixPathTimetables, outputFolder):

    allExistingStopIds = getAllStopNumbers(pathToStopsFile)


    numberOfStops = len(allExistingStopIds)
    counter = 0.0
    dictOfAllStopsData = dict()
    # for stopNumber in sorted(allExistingStopIds):
    for stopNumber in sorted([36234384, 36234327]):
        print( str( 100.0*(counter/numberOfStops)) + "% done (" + str(int(counter)) + " of " + str(numberOfStops) + ")")
        currentStopData = processASingleTimetableForBusStop(stopNumber, prefixPathTimetables, suffixPathTimetables)
        dictOfAllStopsData[stopNumber] = currentStopData
        counter = counter + 1.0

    pickle.dump( dictOfAllStopsData, open( outputFolder + "/timetables.p", "wb" ) )




if __name__ == "__main__":

    pathToFileFromAPIStops = "/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/stops.txt"

    pathPrefixForAllTimetables = "/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/timetables_new/"
    pathSuffixForAllTimetables = ".txt"

    pathToFolderToSavePickledResult = "/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/timetables_pickled_new"





    main(pathToFileFromAPIStops, pathPrefixForAllTimetables, pathSuffixForAllTimetables, pathToFolderToSavePickledResult)




