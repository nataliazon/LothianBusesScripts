import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import operator
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import scipy.stats as stats
import calendar
from pprint import pprint
import datetime
import json
import pickle
import simplejson


import random

set_of_existing_names  =set()


# the below functions needed for reorganizeDataCollected()
def d3l():
    return defaultdict(d2l)


def d2l():
    return defaultdict(dlist)


def d3():
    return defaultdict(d2)


def d2():
    return defaultdict(d1)


def d1():
    return defaultdict(dict)


def dlist():
    return defaultdict(list)

#heatmap exapmle
#
# a = np.random.random((16, 16))
#
# b = [[0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
#      [5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
#      [0, 1, 2, 2, 4, 0, 1, 2, 3, 4],
#      [0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
#      [0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
#      [5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
#      [0, 1, 2, 2, 4, 0, 1, 2, 3, 4],
#      [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
#      ]
#
# plt.imshow(b, cmap='jet', interpolation='nearest')
# plt.show()


# We want to generate heatmap plots at snapshots of time
# For each hour of the day:
# x - bus stop numbers
# y - bus stop numbers
# z - journey times
# so for each snapshot, we have a triangle

def testPlotting():
    a = np.random.random((16, 16))

    exampleData = list()
    for busIdFrom in range (0,126, 1):
        exampleData.append(list())
        for busIDTo in range(0,126,1):
            if busIDTo == busIdFrom:
                val = 0
            elif busIDTo < busIdFrom:
                val = -10
            else:
                val = random.choice([15,28,33,45,12,13,10,1,2,6,5])
            exampleData[busIdFrom].append(val)

    b = [[0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
         [5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
         [0, 1, 2, 2, 4, 0, 1, 2, 3, 4],
         [0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
         [0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
         [5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
         [0, 1, 2, 2, 4, 0, 1, 2, 3, 4],
         [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
         ]

    c = [[0,1,2],
         [0,0,6 ],
         [0, 0, 0]
         ]


    fig, ax = plt.subplots()
    cax = ax.imshow(exampleData, cmap='magma', interpolation='nearest')
    ax.set_title("Journey times")
    cbar = fig.colorbar(cax, ticks=[-10,0,15,30,60])
    cbar.ax.set_yticklabels(["No journey", "same stop", "15 min", "30 min", "1 hour"])

   # plt.imshow(exampleData, cmap='viridis', interpolation='nearest')
   # plt.legend()

    plt.show()


def reorganize_live_data_collected(path_to_data, dataset_name):
    month_number = dict((v, k) for k, v in enumerate(calendar.month_abbr))
    files = os.listdir(path_to_data)
    valid_filename_regexp = re.compile('dataset2_[0-9]{4}-[A-z]{3}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}\.txt')

    valid_files = []

    for filename in files:
        if valid_filename_regexp.match(filename):
            year = int(filename[-24:-20])
            month = int(month_number[filename[-19:-16]])
            day = int(filename[-15:-13])
            hour = int(filename[-12:-10])
            minute = int(filename[-9:-7])
            second = int(filename[-6:-4])
            datetime_stamp = datetime.datetime(year, month, day, hour, minute, second)

            valid_files.append([path_to_data + '/' + filename, datetime_stamp])
    pprint(valid_files)

    reorganized = defaultdict(d3)

    print('Processing ' + str(len(valid_files)) + ' valid data files...')
    for datafile in sorted(valid_files):
        print(datafile)
        data = json.load(open(datafile[0]))
        time_of_request = datafile[1]

        bus_info_data = data['vehicles']
        for busSnapshot in bus_info_data:
            bus_destination = busSnapshot['destination']
            bus_journey_id = busSnapshot['journey_id']
            bus_latitude = busSnapshot['latitude']
            bus_longitude = busSnapshot['longitude']
            bus_next_stop_id = busSnapshot['next_stop_id']
            bus_service_name = busSnapshot['service_name']
            if bus_service_name is None:
                continue
            bus_speed = busSnapshot['speed']
            bus_vehicle_id = busSnapshot['vehicle_id']
            reorganized[bus_service_name][bus_destination][bus_journey_id][time_of_request] = \
                {'latitide': bus_latitude, 'longitude': bus_longitude, 'next_stop_id': bus_next_stop_id,
                 'speed': bus_speed, 'vehicle_id': bus_vehicle_id}

    for serviceName in sorted(reorganized.keys()):
        print(serviceName)
        try:
            pickle_out = open(
                "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/"
                "batch7twoDays/"+dataset_name+"_service" + "_" + serviceName + ".pickle",
                "wb")
            pickle.dump(reorganized[serviceName], pickle_out)
            pickle_out.close()
        except TypeError:
            print("Type error for service ")
            continue

    return

def get_real_time_data(path_to_pickled_real_time_data, service, destination, dataset, timetable_start_in_min, first_stop_ids, second_stop_ids, last_stop_ids, wrong_stop_ids):
    print("Running get_real_time_data ...")
    print("Loading pickle file for service " + service + " to " + destination + "  from path: " + path_to_pickled_real_time_data)
    filename = dataset + "_service_" + service + ".pickle"
    path_to_open_file = path_to_pickled_real_time_data + "/" + filename
    #print(path_to_open_file)
    opened_file = open(path_to_open_file, 'rb')

    loaded_service_pickle = pickle.load(opened_file)
    for item in loaded_service_pickle.keys():
        #print(item)
        pass
    if destination not in loaded_service_pickle.keys():
        raise Exception("Exception. There is no such destination (" + destination + ") for service " + service + ".")

    requested_data = loaded_service_pickle[destination]

    # print(requested_data['6965'].keys())

    # divide by days
    # divide by kind of days (weekdays, weekends, etc)
    # for each journey_ID find first and last seen datetimestamps
    # division: kinds of days / days / journeys

    # datetime.weekday() returns the day of the week
    # 0 = Monday, 6 = Sunday
    # in lothian api day is eithr 0, 6 or 7 - TODO check what do they mean (probably 0-weekdays,6-sat,7-sun)

    weekdays = {0, 1, 2, 3, 4}  # 0
    saturday = {5}  # 6
    sunday = {6}  # 7

    data_by_daytype = defaultdict(list)

    for journey_id in requested_data.keys():
        for datetime_stamp in sorted(requested_data[journey_id]):
            #print("For the date " + str(datetime_stamp) + " the weekday is " + str(datetime_stamp.weekday()))
            if datetime_stamp.weekday() in weekdays:
                data_by_daytype[0].append([journey_id, datetime_stamp, requested_data[journey_id][datetime_stamp]])
            elif datetime_stamp.weekday() in saturday:
                data_by_daytype[0].append([journey_id, datetime_stamp, requested_data[journey_id][datetime_stamp]])
            elif datetime_stamp.weekday() in sunday:
                data_by_daytype[0].append([journey_id, datetime_stamp, requested_data[journey_id][datetime_stamp]])
            else:
                raise Exception("Invalid weekday returned for timestamp")
    # print(data_by_daytype.keys())

    data_by_daytype_by_day = defaultdict(dlist)

    for type_of_day in data_by_daytype.keys():
        for datetime_stamp in data_by_daytype[type_of_day]:
            # print(datetime_stamp)
            day_number = (datetime_stamp[1]).date().day
            month_number = (datetime_stamp[1]).date().month
            year_number = (datetime_stamp[1]).date().year
            date_value = datetime.date(year_number, month_number, day_number)
            data_by_daytype_by_day[type_of_day][date_value].append(datetime_stamp)

    # only process full days
    # In order to do that, find first and last timeStamp for each day

    acceptable_earliest_time = datetime.time(8, 0, 0, 0)
    acceptable_latest_time = datetime.time(20, 59, 59, 59)

    # days having timestamps at least between acceptable_earliest and acceptable_latest time values
    full_day_dates = set()

    for type_of_day in sorted(data_by_daytype_by_day.keys()):
        for date_value in sorted(data_by_daytype_by_day[type_of_day]):
            earliest_time = datetime.time(23, 59, 59, 59)
            latest_time = datetime.time(0, 0, 0, 0)
            print("Processing day " + str(date_value))

            for snapshot_data in data_by_daytype_by_day[type_of_day][date_value]:
                #print(snapshot_data[1])
                if snapshot_data[1].time() < earliest_time:
                    earliest_time = snapshot_data[1].time()
                if snapshot_data[1].time() > latest_time:
                    latest_time = snapshot_data[1].time()
            #print("Earliest " + str(earliest_time))
            #print("Latest " + str(latest_time))
            if earliest_time < acceptable_earliest_time and latest_time > acceptable_latest_time:
                full_day_dates.add(date_value)

    print("Full days: " + str(full_day_dates))
    data_by_daytype_by_day_by_journey = defaultdict(d2l)

    for type_of_day in sorted(data_by_daytype_by_day.keys()):
        for date_value in sorted(data_by_daytype_by_day[type_of_day]):
            if date_value in full_day_dates:
                for snapshot_data in data_by_daytype_by_day[type_of_day][date_value]:
                    data_by_daytype_by_day_by_journey[type_of_day][date_value][snapshot_data[0]].append(
                        [snapshot_data[1], snapshot_data[2]])

    #data_by_daytype_by_day_by_journey


    #return data_by_daytype_by_day_by_journey

    new_data_by_daytype_by_day_by_journey = defaultdict(d2l)


    first_and_last_seen = dict()
    for type_of_day in sorted(data_by_daytype_by_day_by_journey.keys()):
        first_and_last_seen[type_of_day] = dict()
        for date_value in sorted(data_by_daytype_by_day_by_journey[type_of_day].keys()):
            first_and_last_seen[type_of_day][date_value] = dict()
            for journey_id in sorted(data_by_daytype_by_day_by_journey[type_of_day][date_value].keys()):
                first_and_last_seen[type_of_day][date_value][journey_id] = dict()
                #print(journey_id)
                # print(data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id])

                first_seen_snapshot = min(data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id],
                        key=operator.itemgetter(0))
                last_seen_snapshot = max(data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id],
                        key=operator.itemgetter(0))

                #print("F_S: " + str(first_seen_snapshot))
                #print("L_S: " + str(last_seen_snapshot))

                if(first_seen_snapshot[1]['next_stop_id'] in first_stop_ids or first_seen_snapshot[1]['next_stop_id']in second_stop_ids):
                    new_data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id] = data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id]
                elif (first_seen_snapshot[1]['next_stop_id'] in wrong_stop_ids or first_seen_snapshot[1]['next_stop_id']in wrong_stop_ids or last_seen_snapshot[1]['next_stop_id'] in wrong_stop_ids):
                    print("Nonexistant stop bus " + str(journey_id) + " first " + str(
                       first_seen_snapshot[1]['next_stop_id']) + " last " + last_seen_snapshot[1]['next_stop_id'])
                    continue
                else:
                    print("Misbehaved bus " + str(journey_id) + " first " + str(first_seen_snapshot[1]['next_stop_id']) + " last " + last_seen_snapshot[1]['next_stop_id'])
                    continue

    return new_data_by_daytype_by_day_by_journey




#from services we need:
#for each stop id obtain the previous stop id
#for each stop if obtain it's number in the sequence that a journey follow

def get_stops_data(path_to_stops_file):
    fullFile = open(path_to_stops_file).read()

    beg = fullFile.find("<")
    end = fullFile.rfind(">")

    correctJsonFile = fullFile[end + 1:]

    stopsLoaded = simplejson.loads(correctJsonFile)

    stopsById = dict()

    for item in stopsLoaded['stops']:
        #print(item)
        stopsById[item['stop_id']] = item

    return stopsById

def get_services_routes_data(path_to_services_file, stop_data_by_stop_id):


    fullFile = open(path_to_services_file).read()

    beg = fullFile.find("<")
    end = fullFile.rfind(">")

    correctJsonFile = fullFile[end + 1:]

    servicesLoaded = simplejson.loads(correctJsonFile)


    list_all_services = servicesLoaded['services']

    bus_stops_by_dest_by_service = dict()
    for service in list_all_services:
        #print(service)
        routes = service['routes']
        bus_stops_by_routes = dict()
        for route in routes:
            points = route['points']
            if service['name'] == "5":
                pass
                #print(service['name'])
                #print(route['destination'])
                #print(points)
            #print("-------------")
            bus_stop_counter = 0
            bus_stopped_mapped_to_indexes = dict()
            for routePoint in points:
                if routePoint['stop_id'] != None and routePoint['stop_id'] in stop_data_by_stop_id.keys():
                    #print(routePoint)
                    bus_stopped_mapped_to_indexes[bus_stop_counter] = [True,routePoint['stop_id'], [routePoint['latitude'], routePoint['longitude']],stop_data_by_stop_id[routePoint['stop_id']]['name']]
                    bus_stop_counter +=1
                elif routePoint['stop_id'] != None:
                    #bus_stopped_mapped_to_indexes[bus_stop_counter] = ["?", routePoint['stop_id'],[routePoint['latitude'], routePoint['longitude']],names_generator(1)]
                    print("!! Service " + str(service['name']) + " to " + str(route['destination']) + " has nonexistant stop: " + str(routePoint['stop_id']) + " skipping ...")
                    #bus_stop_counter += 1
                else:
                    bus_stopped_mapped_to_indexes[bus_stop_counter] = [False, bus_stop_counter+42000000,[routePoint['latitude'], routePoint['longitude']], names_generator(1)]
                    bus_stop_counter+=1
            bus_stops_by_routes[route['destination']] = bus_stopped_mapped_to_indexes
        bus_stops_by_dest_by_service[service['name']] = bus_stops_by_routes

    # print(bus_stops_by_dest_by_service)

    return bus_stops_by_dest_by_service


def process_single_journey(single_journy_data):
    print(single_journy_data)
    data_by_next_stop_id=defaultdict(list)
    for item in single_journy_data:
        data_by_next_stop_id[item[1]['next_stop_id']].append(item)
    print(data_by_next_stop_id.keys())
    for next_stop_id in data_by_next_stop_id.keys():
        print(next_stop_id)
        print("snapshots: " + str(len(data_by_next_stop_id[next_stop_id])))
        print(data_by_next_stop_id[next_stop_id])


    return


def get_rawjavadata(path_to_txt_file_with_data):
    print("Running get_rawjavadata ...")
    opened_datafile = open(path_to_txt_file_with_data)
    linesofdata = opened_datafile.readlines()
    #

    # read only the appropriate data values: busID, time, start_time
    datalist = []
    for lineofdata in linesofdata:
        # print(lineofdata)
        splat = lineofdata.split()
        # splat[2]==bID , splat[4]==currentNow, splat[6]==bus_start_time
        datalist.append([int(splat[2]), float(splat[4]), float(splat[6])])

    # This data contains multiple simulation instances (replications).
    # The code below divides the data into sets depending on the replication ID

    results_per_simuation_instance = defaultdict(list)
    last_timeof_simulation = 0.0
    current_instance = 0
    for sublist in datalist:
        # print(sublist)
        time_of_simulation = sublist[1]
        if time_of_simulation < last_timeof_simulation:
            current_instance = current_instance + 1

        results_per_simuation_instance[current_instance].append(sublist)
        last_timeof_simulation = time_of_simulation

    # For testing purposes:
    # for instanceNumber in results_per_simulation_instance.keys():
    # print(instanceNumber)
    # for item in results_per_simulation_instance[instanceNumber]:
    #    print(str(instanceNumber) +  "    "+str(item[1]))

    # now, for each simulation instance we need to count buses that are active at every timestep of the simulation
    # 1. For each busID find lastSeen time (max(timeOfSnap))
    # 2. Assume the bus was active between busStart and lastSeen time
    # 3. For each timevalue count buses for whom busStart < timevalue < lastSeen

    results_per_simulation_instance_per_bus_id = dict()
    for simulation_instance in results_per_simuation_instance.keys():
        this_instance_results = results_per_simuation_instance[simulation_instance]
        snapshots_per_bus_id = defaultdict(list)
        for snapshot in this_instance_results:
            snapshots_per_bus_id[snapshot[0]].append(snapshot)
        results_per_simulation_instance_per_bus_id[simulation_instance] = snapshots_per_bus_id

    # results_per_simulation_instance_per_bus_i_d[simulation_instance][busID] = [[bus_id, timeOfSnap, start_time], .. ]

    bus_life_ranges_per_simulation_id_per_bus_id = defaultdict(dict)

    for simulationIns in results_per_simulation_instance_per_bus_id.keys():
        for bus_id in results_per_simulation_instance_per_bus_id[simulationIns].keys():
            # print(results_per_simulation_instance_per_bus_i_d[simulationIns][bus_id])

            last_seen_snapshot = (
                max(results_per_simulation_instance_per_bus_id[simulationIns][bus_id], key=operator.itemgetter(1)))
            bus_start_time = last_seen_snapshot[2]
            bus_last_seen_time = last_seen_snapshot[1]
            # print("instance: " + str(simulationIns) + " busID: " + str(bus_id) + " start_time: " + str(bus_start_time)
            # + " last_seen_time: " + str(bus_last_seen_time))
            bus_life_ranges_per_simulation_id_per_bus_id[simulationIns][bus_id] = [bus_start_time, bus_last_seen_time]

    # For each simulation instance, count the buses active at every simulation step.

    active_and_inactive_buses_per_simulation_id_per_time = defaultdict(dict)

    for simulation_instance_number in bus_life_ranges_per_simulation_id_per_bus_id.keys():
        for i in range(0, 4500, 1):
            current_simulation_time = float(i)
            number_of_active_buses = 0
            number_of_inactive_buses = 0
            for bus_id in bus_life_ranges_per_simulation_id_per_bus_id[simulation_instance_number].keys():
                this_bus_data = bus_life_ranges_per_simulation_id_per_bus_id[simulation_instance_number][bus_id]
                start_time = this_bus_data[0]
                last_seen_time = this_bus_data[1]
                if start_time <= current_simulation_time < last_seen_time:
                    number_of_active_buses = number_of_active_buses + 1
                else:
                    number_of_inactive_buses = number_of_inactive_buses + 1
            active_and_inactive_buses_per_simulation_id_per_time[simulation_instance_number][
                current_simulation_time] = [
                number_of_active_buses, number_of_inactive_buses]

    # print(active_and_inactive_buses_per_simulation_id_per_time)

    simulation_times = []

    for i in range(0, 4500, 1):
        current_simulation_time = int(i)
        simulation_times.append(current_simulation_time)

    active_buses_values_list_per_times = defaultdict(list)
    inactive_buses_values_list_per_times = defaultdict(list)

    for simulation_instance in active_and_inactive_buses_per_simulation_id_per_time.keys():
        for simulationTime in simulation_times:
            active_b = active_and_inactive_buses_per_simulation_id_per_time[simulation_instance][simulationTime][0]
            inactive_b = active_and_inactive_buses_per_simulation_id_per_time[simulation_instance][simulationTime][1]
            active_buses_values_list_per_times[simulationTime].append(active_b)
            inactive_buses_values_list_per_times[simulationTime].append(inactive_b)
            # print(str(simulation_instance)+", " + str(simulationTime) + ", " + str(active_b) + ", " + str(inactive_b))
            pass

    # compute the mean, standard deviation and  standard error of the mean and equal to std / sqrt(n) where n is the
    # sample size. lets assume that the value of delta degrees of freedom = 1\

    # numpy std has the default value of deltaDegreesOfFreedom ddof=0, and scipy.sem has the default value ddof=1
    # which should we use? ?? ?? -> to be established, right now using 0  because we want to divide by sqrt(N-0),
    # so that our sem is std/sqrt(N) this just means that there is no additional constraints on the data (?)

    # Functions for computing the mean, std and sem print("numpy mean " + str(np.mean([1,2,3,4,5,6,7,8,9,
    # 10]))) print("numpy std " + str(np.std([1,2,3,4,5,6,7,8,9,10]))) print("scipy.stats sem " + str(stats.sem([1,2,
    # 3,4,5,6,7,8,9,10], ddof=1))) #standard error of the mean print("manual sem " + str(np.std([1,2,3,4,5,6,7,8,9,
    # 10], ddof=1)/np.sqrt(len([1,2,3,4,5,6,7,8,9,10])))) #standard error of the mean
    #

    active_buses_mean_std_sem_per_simulation_time = dict()

    for simulationTime in active_buses_values_list_per_times.keys():
        list_of_buses_counts = active_buses_values_list_per_times[simulationTime]
        # print(list_of_buses_counts)
        mean = np.mean(list_of_buses_counts)
        std = np.std(list_of_buses_counts, ddof=0)
        sem = stats.sem(list_of_buses_counts, ddof=0)
        active_buses_mean_std_sem_per_simulation_time[simulationTime] = [mean, std, sem]

    inactive_buses_mean_std_sem_per_simulation_time = dict()

    for simulationTime in inactive_buses_mean_std_sem_per_simulation_time.keys():
        list_of_buses_counts = inactive_buses_mean_std_sem_per_simulation_time[simulationTime]
        # print(list_of_buses_counts)
        mean = np.mean(list_of_buses_counts)
        std = np.std(list_of_buses_counts, ddof=0)
        sem = stats.sem(list_of_buses_counts, ddof=0)
        inactive_buses_mean_std_sem_per_simulation_time[simulationTime] = [mean, std, sem]

    list_of_simulation_times = list()
    list_of_active_means = list()
    list_of_active_stds = list()
    list_of_active_sems = list()
    list_of_inactive_means = list()
    list_of_inactive_stds = list()
    list_of_inactive_sems = list()

    for simulationTime in active_buses_mean_std_sem_per_simulation_time.keys():
        list_of_simulation_times.append(simulationTime)
        list_of_active_means.append(active_buses_mean_std_sem_per_simulation_time[simulationTime][0])
        list_of_active_stds.append(active_buses_mean_std_sem_per_simulation_time[simulationTime][1])
        list_of_active_sems.append(active_buses_mean_std_sem_per_simulation_time[simulationTime][2])

    for simulationTime in inactive_buses_mean_std_sem_per_simulation_time.keys():
        list_of_inactive_means.append(inactive_buses_mean_std_sem_per_simulation_time[simulationTime][0])
        list_of_inactive_stds.append(inactive_buses_mean_std_sem_per_simulation_time[simulationTime][1])
        list_of_inactive_sems.append(inactive_buses_mean_std_sem_per_simulation_time[simulationTime][2])

    print("Returning from get_rawjavadata ... ")
    return {"simulation_times": list_of_simulation_times, "activeBusesMeans": list_of_active_means,
            "activeBusesStds": list_of_active_stds, "activeBusesSems": list_of_active_sems,
            "inactiveBusesMeans": list_of_inactive_means, "inactiveBusesStds": list_of_inactive_stds,
            "inactiveBusesSems": list_of_inactive_sems}


def get_simulation_data(path_to_data):
    print("Getting simulation data...." + path_to_data)

    opened_datafile = open(path_to_data)
    linesofdata = opened_datafile.readlines()
    #

    # read only the appropriate data values: busID, time, start_time
    datalist = []
    for lineofdata in linesofdata:
        #print(lineofdata)
        splat = lineofdata.split()
        # splat[2]==bID , splat[4]==currentNow, splat[6]==bus_start_time

        date = str(splat[0])
        bId = int(splat[2])
        currentNow = float(splat[4])
        busStartTime = float(splat[6])
        busLocationId = int(splat[8])
        busX = float(splat[10])
        busY = float(splat[12])

        datalist.append({'datetime':date, 'busID':bId, 'now':currentNow, 'startTime':busStartTime, 'loc':busLocationId, 'x':busX, 'y':busY})
        #datalist.append([int(splat[2]), float(splat[4]), float(splat[6])])

    # This data contains multiple simulation instances (replications).
    # The code below divides the data into sets depending on the replication ID

    results_per_simuation_instance = defaultdict(list)
    last_timeof_simulation = 0.0
    current_instance = 0
    for sublist in datalist:
        #print(sublist)

        #for range, stopfrom/stopTo, list of journeytimes
        #
        time_of_simulation = sublist['now']
        if time_of_simulation < last_timeof_simulation:
            current_instance = current_instance + 1

        results_per_simuation_instance[current_instance].append(sublist)
        last_timeof_simulation = time_of_simulation

    #print(results_per_simuation_instance.keys())
    #print(results_per_simuation_instance[0][0])


    return results_per_simuation_instance


def get_stop_number_from_carma_id(location):
    # TODO

    pass


def process_simulation_data(results_per_simuation_instance, pairs, pickle_out):
    print("Running process simulation data ...")
    print(pairs)
    print(results_per_simuation_instance.keys())

    ranges = define_ranges()
    hour_ranges = ranges[0]
    min_ranges = ranges[1]

    #hour_ranges= [[13.5, 14.5]]

    data_per_range = dict()
    for current_range in min_ranges:
        print("Processing " + str(current_range))
        range_as_str = str(current_range[0]) + "-" + str(current_range[1])
        time_from = current_range[0]
        time_to = current_range[1]
        data_per_range[range_as_str] = dict()
        for fromId in pairs.keys():
            data_per_range[range_as_str][fromId] = dict()
            for toId in pairs[fromId].keys():
                print("   Processing " + str(current_range) + " / " + str(fromId) + " to " + str(toId))

                data_per_range[range_as_str][fromId][toId]=dict()
                dataFrom = pairs[fromId][toId][0]
                dataTo = pairs[fromId][toId][1]
                #print(fromId)
                from_Carma_id = dataFrom[0][1]
                #print("*" + str(from_Carma_id))
                to_Carma_id = dataTo[0][1]
                #print(toId)
                #print("&" +str(to_Carma_id))
                matching_events_per_sim_instance= dict()
                journey_times = list()
                #print("number of results  " + str(len(results_per_simuation_instance.keys())))
                for simulation_instance in results_per_simuation_instance.keys():
                    list_of_collected_data_for_current_instance = results_per_simuation_instance[simulation_instance]
                    list_of_collected_data_for_current_instance_per_bus_id = defaultdict(list)
                    for item in list_of_collected_data_for_current_instance:
                        bus_id = item['busID']
                        list_of_collected_data_for_current_instance_per_bus_id[bus_id].append(item)
                    for item in list_of_collected_data_for_current_instance:
                        #print(item)
                        location = item['loc']
                        location_translated = get_stop_number_from_carma_id(location)
                        time = item['now']
                        bus_id = item['busID']
                        if(time_from < time  < time_to):
                            #print("time_from < time  < time_to")
                            #print("loc " + str(location) + "   toid" + str(toId) +  "   fromid" + str(fromId))
                            if(location == toId):
                                #print("time_from < time  < time_to and location == toId")
                                all_snapshots_of_this_bus = list_of_collected_data_for_current_instance_per_bus_id[bus_id]
                                for snap in all_snapshots_of_this_bus:
                                    if snap['loc'] == fromId:

                                        journey_time = time - snap['now']
                                        journey_times.append(journey_time)
                average_journey_time = None
                if len(journey_times) > 0:
                    print("Journey times len " + str(len(journey_times)))

                    average_journey_time = np.average(journey_times)
                    print("avg: " + str(average_journey_time))
                print("Average journey time :  " + str(average_journey_time))
                data_per_range[range_as_str][fromId][toId] = average_journey_time


    try:
        pickle_out = open(
            pickle_out,
            "wb")
        pickle.dump(data_per_range, pickle_out)
        pickle_out.close()
    except TypeError:
        print("Type error for service ")



    return data_per_range


def get_all_bus_stop_pairs(bus_stops_by_route_by_service, service, destination, include_non_stop_points=True):

    #print(bus_stops_by_route_by_service.keys())
    #print(bus_stops_by_route_by_service['5'].keys())
    #print(bus_stops_by_route_by_service[service][destination])

    bus_stop_pairs = dict()


    if include_non_stop_points:
        for stop_number_from in (bus_stops_by_route_by_service[service][destination].keys()):
            bus_stop_pairs[stop_number_from] = dict()
            for stop_number_to in (bus_stops_by_route_by_service[service][destination].keys()):
                if stop_number_to <= stop_number_from:
                   continue
                   # print("From " + str(stop_number_from) + " to " + str(stop_number_to) + " is 0")
                else:
                    print("From " + str(stop_number_from) + " to " + str(stop_number_to))


                    point_from = bus_stops_by_route_by_service[service][destination][stop_number_from]
                    point_to = bus_stops_by_route_by_service[service][destination][stop_number_to]

                    point_from_with_number = [[stop_number_from], point_from]
                    point_to_with_number = [[stop_number_to], point_to]

                    bus_stop_pairs[stop_number_from][stop_number_to] = [point_from_with_number, point_to_with_number]

    else:

        #remove non-stops from data
        new_data_without_non_stops = dict()
        counter = 0
        for old_routepoint_number in bus_stops_by_route_by_service[service][destination].keys():
            current_point = bus_stops_by_route_by_service[service][destination][old_routepoint_number]
            if(current_point[0] == True) and old_routepoint_number != 83:
                if old_routepoint_number<83:
                    new_point = [[counter, old_routepoint_number], current_point]
                else:
                    new_point = [[counter, old_routepoint_number-1], current_point]
                if counter == 39:
                    print("%%%%" + str(bus_stops_by_route_by_service[service][destination][old_routepoint_number]))


                new_data_without_non_stops[counter] = new_point
                counter+=1
            else:
                pass
        #print(new_data_without_non_stops)
        for no_from in new_data_without_non_stops.keys():
            bus_stop_pairs[no_from] = dict()
            for no_to in new_data_without_non_stops.keys():
                if no_to <= no_from:
                   continue
                else:
                    bus_stop_pairs[no_from][no_to] = [new_data_without_non_stops[no_from], new_data_without_non_stops[no_to]]

    #print(bus_stop_pairs)


    return bus_stop_pairs


def define_ranges():
    start = 0.5
    finish = 23.5

    #start=17.5
    #finish = 18.5

    hour_ranges = []
    min_ranges = []
    while start < finish:

        hour_ranges.append([start, start+1])
        min_ranges.append([start*60.0, (start+1)*60.0])
        start += 1

    return [hour_ranges, min_ranges]


def process_real_time_data(real_time_data, all_pairs, pickle_out):

    ranges = define_ranges()

    hour_ranges = ranges[0]
    min_ranges = ranges[1]

    print(hour_ranges)
    print(len(hour_ranges))
    print(min_ranges)

    #for each day
    # for each pair
    # for each hour period

    counter_of_instances = 0
    times_per_daytype_day_journey_range_from_to = dict()
    for day_type in real_time_data.keys():
        current_day_journey_times = dict()
        for day in real_time_data[day_type].keys():
            current_journey_journey_times = dict()
            for journey in real_time_data[day_type][day].keys():
                current_range_journey_times = dict()
                for time_range in min_ranges:
                    current_pair_journey_times = dict()
                    hash_range = str(time_range[0]) + "-" + str(time_range[1])
                    for pair_from in all_pairs.keys():
                        current_pair_journey_times[pair_from] = dict()
                        for pair_to in all_pairs[pair_from].keys():
                            current_pair_journey_times[pair_from][pair_to] = None
                            current_journey_data = real_time_data[day_type][day][journey]
                            if len(current_journey_data) == 0:
                                print("Len of current_journey_data is 0 for " + str(pair_from) +  "  " + str(pair_to))
                                continue
                            elif len(current_journey_data[0]) ==0:
                                print("Len of current_journey_data[0] is 0 for " + str(pair_from) + "  " + str(pair_to))
                                continue
                            current_journey_start = current_journey_data[0][0]
                            current_journey_finish = current_journey_data[len(current_journey_data)-1][0]
                            start_min = current_journey_start.minute + current_journey_start.hour*60.0
                            fin_min = current_journey_finish.minute + current_journey_finish.hour*60.0
                            counter_of_instances += 1
                            if fin_min < time_range[0] or start_min > time_range[1]:
                                continue
                            else:
                                stop_id_to = all_pairs[pair_from][pair_to][1][1][1]
                                stop_id_from = all_pairs[pair_from][pair_to][0][1][1]
                                #print("Stop to " + str(pair_to) + "(" + str(stop_id_to) + ")")
                                snaps_with_stop_id_to_as_next_stop = []
                                snaps_with_stop_id_from_as_next_stop = []
                                for timesnap in current_journey_data:
                                    timesnap_next_stop = timesnap[1]['next_stop_id']
                                    if(str(timesnap_next_stop) == str(stop_id_to)):
                                        snaps_with_stop_id_to_as_next_stop.append(timesnap)
                                    if(str(timesnap_next_stop) == str(stop_id_from)):
                                        snaps_with_stop_id_from_as_next_stop.append(timesnap)
                                #print("This many occurences found to:" + str(len(snaps_with_stop_id_to_as_next_stop)))
                                #print("This many occurences found from:" + str(len(snaps_with_stop_id_from_as_next_stop)))
                                if(len(snaps_with_stop_id_to_as_next_stop) == 0 or len(snaps_with_stop_id_from_as_next_stop) == 0):
                                    continue
                                arrival_at_stop_to = snaps_with_stop_id_to_as_next_stop[-1]
                                arrival_at_stop_from = snaps_with_stop_id_from_as_next_stop[-1]
                                #example [datetime.datetime(2018, 3, 6, 7, 18, 1), {'latitide': 55.955417, 'longitude': -3.154888, 'next_stop_id': '36234737', 'speed': 16, 'vehicle_id': '951'}]
                                minutes_of_arrival_at_stop_to = arrival_at_stop_to[0].hour*60.0 + arrival_at_stop_to[0].minute + arrival_at_stop_to[0].second/60.0
                                minutes_of_arrival_at_stop_from = arrival_at_stop_from[0].hour*60.0 + arrival_at_stop_from[0].minute + arrival_at_stop_from[0].second/60.0
                                journey_time = minutes_of_arrival_at_stop_to - minutes_of_arrival_at_stop_from
                                #print("JOURNEY TIME: " + str(journey_time))
                                current_pair_journey_times[pair_from][pair_to] = journey_time
                                #print(arrival_at_stop_to)
                                #print(minutes_of_arrival_at_stop_to)
                                print(str(100 * counter_of_instances / 13694844) + " %  done")
                    current_range_journey_times[hash_range] = current_pair_journey_times




                                    #print("Zero occurences of next stop " + str(pair_to) + "(" + str(stop_id_to) + ")" + "for journey " + str(current_journey_data))
                                    #for item in current_journey_data:
                                     #   print(item[1]['next_stop_id'])



                current_journey_journey_times[journey] = current_range_journey_times
            current_day_journey_times[day] = current_journey_journey_times
        times_per_daytype_day_journey_range_from_to[day_type] = current_day_journey_times

    #print(times_per_daytype_day_journey_range_from_to)

    try:
        pickle_out = open(
            pickle_out,
            "wb")
        pickle.dump(times_per_daytype_day_journey_range_from_to, pickle_out)
        pickle_out.close()
    except TypeError:
        print("Type error for service ")

    # find journeys that are during that hour at the to_stop
    # get their time the from_stop
    # find the difference => journey time
    # put into journey_length_per_type_of_day_per_day_per_hour_from_to = X
    return


def read_pickled_data_realtime(real_time_data_pickle_path):
    opened_file = open(real_time_data_pickle_path, 'rb')

    loaded_pickle = pickle.load(opened_file)

    #print(loaded_pickle[0][datetime.date(2018, 3, 7)]['6804']['750.0-810.0'][12][15])

    data_by_range = dict()

    for day in loaded_pickle[0].keys():
        for journey in loaded_pickle[0][day].keys():
            for range in loaded_pickle[0][day][journey].keys():
                if range not in data_by_range.keys():
                    data_by_range[range] = dict()
                for from_stop in loaded_pickle[0][day][journey][range].keys():
                    if from_stop not in data_by_range[range].keys():
                        data_by_range[range][from_stop] = dict()
                    for to_stop in loaded_pickle[0][day][journey][range][from_stop].keys():
                        if to_stop not in data_by_range[range][from_stop].keys():
                            data_by_range[range][from_stop][to_stop] = list()
                        if loaded_pickle[0][day][journey][range][from_stop][to_stop] != None:
                            data_by_range[range][from_stop][to_stop].append(loaded_pickle[0][day][journey][range][from_stop][to_stop])
                        else:
                            if(to_stop == 39 or from_stop == 39) :
                                pass
                                #print("None data for: " + str(day) + "/" + str(journey) + "/" +str(range) + "/" + str(from_stop) + "/" + str(to_stop))

    print(str(data_by_range))

    average_journey_time_by_range_by_from_by_to = dict()

    for range in data_by_range.keys():
        average_journey_time_by_range_by_from_by_to[range] = dict()
        for from_stop in data_by_range[range].keys():
            average_journey_time_by_range_by_from_by_to[range][from_stop] = dict()
            for to_stop in data_by_range[range][from_stop].keys():
                list_of_journey_times = data_by_range[range][from_stop][to_stop]
                if(len(list_of_journey_times)) > 0:
                    average_value = np.mean(list_of_journey_times)
                else:
                    average_value = None
                average_journey_time_by_range_by_from_by_to[range][from_stop][to_stop] = average_value


    return average_journey_time_by_range_by_from_by_to


def plot_diff(differences, real_data_vals, sim_vals, path_to_savefig, title=""):
    print(differences)

    a = list()
    b= list()
    c=list()

    for stop_from in differences.keys():
        a.insert(stop_from,list())
        b.insert(stop_from,list())
        c.insert(stop_from,list())
        for stop_to in differences.keys():
            journey_time = None
            real_jt = None
            sim_jt = None
            if stop_to in differences[stop_from].keys():
                journey_time = differences[stop_from][stop_to]
                real_jt = real_data_vals[stop_from][stop_to]
                sim_jt = sim_vals[stop_from][stop_to]
            else:
                journey_time = -3.0

            #print(str(stop_from) + "   " + str(stop_to) + "  " + str(journey_time))

            a[stop_from].insert(stop_to, journey_time)
            b[stop_from].insert(stop_to, real_jt)
            c[stop_from].insert(stop_to, sim_jt)

    b = list()

    fig, ax = plt.subplots()
    counter = 0
    for item in a:
        #print("row:" + str(len(item)))
        b.append(list())

        subcounter = 0
        for subitem in item:
            newSubitem = subitem
            if subitem == None:
                newSubitem = -3.0
                #floatized = float(subitem)
            b[counter].append(newSubitem)
            subcounter+=1
        counter += 1
    #print("column " + str(len(a)))

    middle = title.find("-")
    from_title = str(float(title[0:middle])/60.0)
    to_title = str(float(title[middle+1:])/60.0)



    cax = ax.imshow(b, cmap='magma', interpolation='nearest')
    ax.set_title("Journey times diff")

    cbar = fig.colorbar(cax, ticks=[-3, 1,2,3,4,  5, 7,  15, 20, 25, 30, 45, 60, 360, 1440])
    cbar.ax.set_yticklabels(["No data","1 min", "2 min", "3 min","4 min", "5 min", "7 min", "15 min","20 min", "25 min", "30 min","45 min", "1 hour", "3 hours", "1 day"])
    plt.title(from_title + "  to  " + to_title)



    #plt.show()
    plt.savefig(path_to_savefig)
    plt.close()
    return

def plot_real_data(average_journey_times_real_data, path_to_savefig, title=""):
    print("Processin ")

    print(average_journey_times_real_data)

    a = list()

    for stop_from in average_journey_times_real_data.keys():
        a.insert(stop_from,list())
        for stop_to in average_journey_times_real_data.keys():
            journey_time = None
            if stop_to in average_journey_times_real_data[stop_from].keys():
                journey_time = average_journey_times_real_data[stop_from][stop_to]
            else:
                journey_time = -3.0

            #print(str(stop_from) + "   " + str(stop_to) + "  " + str(journey_time))

            a[stop_from].insert(stop_to, journey_time)

    b = list()

    fig, ax = plt.subplots()
    counter = 0
    for item in a:
        #print("row:" + str(len(item)))
        b.append(list())

        subcounter = 0
        for subitem in item:
            newSubitem = subitem
            if subitem == None:
                newSubitem = -3.0
                #floatized = float(subitem)
            b[counter].append(newSubitem)
            subcounter+=1
        counter += 1
    #print("column " + str(len(a)))

    middle = title.find("-")
    from_title = str(float(title[0:middle])/60.0)
    to_title = str(float(title[middle+1:])/60.0)



    cax = ax.imshow(b, cmap='nipy_spectral', interpolation='nearest')
    ax.set_title("Journey times")

    cbar = fig.colorbar(cax, ticks=[-3, 1,  5, 7,  15, 20, 25, 30, 45, 60, 360, 1440])
    cbar.ax.set_yticklabels(["No journey","1 min", "5 min", "7 min", "15 min","20 min", "25 min", "30 min","45 min", "1 hour", "3 hours", "1 day"])
    plt.title(from_title + "  to  " + to_title)



    #plt.show()
    plt.savefig(path_to_savefig)
    plt.close()
    return

def get_diff_sum(reald, simd):
    sum = 0.0
    sum_data = 0.0
    sum_sim = 0.0
    counter = 0.0
    for fromid in reald.keys():

        for toid in reald[fromid].keys():
            a = reald[fromid][toid]
            b = simd[fromid][toid]
            if(a != None and b!=None):
                diffval = a - b
                diffval = np.sqrt(diffval * diffval)
                sum += diffval
                sum_data += a
                sum_sim += b
                counter += 1.0
            else:
                diffval = None

    if counter != 0:
        return [sum/counter,sum_data/counter, sum_sim/counter]
    else:
        return -3.0


def difference_between_heatmaps(reald, simd):
    diffdata = dict()
    for fromid in reald.keys():
        diffdata[fromid] = dict()
        for toid in reald[fromid].keys():
            a = reald[fromid][toid]
            b = simd[fromid][toid]
            if(a != None and b!=None):
                diffval = a - b
                diffval = np.sqrt(diffval * diffval)
            else:
                diffval = None
            diffdata[fromid][toid] = diffval
    return diffdata


def plot_diff_per_range(difference_per_range, path, pathToTxt):

    x= []
    y = []
    z = []
    p = []
    counter=0
    for keyval in difference_per_range[0].keys():
        x.append(counter)
        y.append(difference_per_range[0][keyval])
        z.append(difference_per_range[1][keyval])
        p.append(difference_per_range[2][keyval])
        counter+=1

    openedFile = open(pathToTxt, "w")
    for i in range(0, len(x)):
        openedFile.write(str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + " " + str(p[i]) + "\n")
    openedFile.close()

    plt.close()
    plt.plot(x, y)
    plt.savefig(path)
    return

def main():


    current_moverate = "traffic"

    path_to_services_file = "/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/data_downloaded_from_api_8_march_2018/services.txt"
    path_to_stops_file = "/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/data_downloaded_from_api_8_march_2018/stops.txt"


    path_to_simulation_data = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/carmadata/after_removing_39/moverate"+current_moverate+"/data.txt"


    path_to_new_sim_with_traffic = '/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/carmadata/after_removing_39/newTraffic/withTraffic.txt'
    path_to_new_sim_without_traffic = '/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/carmadata/after_removing_39/newTraffic/noTraffic.txt'

    stop_data_by_stop_id = get_stops_data(path_to_stops_file)

    bus_stops_by_route_by_service = get_services_routes_data(
        path_to_services_file, stop_data_by_stop_id)

    #get_real_time_data("/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/batch7twoDays", "5", "Hunters Tryst", "7")


    # The below code allows to process all points on a given route like bus stops = it adds bus_stop_ids and names to those unlabelled ones.
    # service_5_hunters_tryst = bus_stops_by_route_by_service["5"]["Hunter's Tryst"]
    # print(service_5_hunters_tryst)
    # for key_of_stop in service_5_hunters_tryst.keys():
    #     print(str(key_of_stop) + " " + str(service_5_hunters_tryst[key_of_stop]))





    #print(names_generator())

    #For service 5, generate pairs of bus stops

    print("=====")
    service = "5"
    destination = "Hunter's Tryst"
    all_pairs = get_all_bus_stop_pairs(bus_stops_by_route_by_service, service, destination, include_non_stop_points=False)

    print("??" + str(all_pairs))
    #For every journey, generate times at bus stops by their IDs

    #sim_data = get_simulation_data(path_to_simulation_data)

    simTraffic = get_simulation_data(path_to_new_sim_with_traffic)
    simNoTraffic = get_simulation_data(path_to_new_sim_without_traffic)

    sim_pickleTraffic = "/Users/natalia/pickledTrafficSim"+current_moverate+".pickle"
    sim_pickleNoTraffic = "/Users/natalia/pickledNoTrafficSim"+current_moverate+".pickle"

    #process_simulation_data(sim_data, all_pairs, sim_pickle)
    #process_simulation_data(simTraffic, all_pairs, sim_pickleTraffic)
    #process_simulation_data(simNoTraffic, all_pairs, simNoTraffic)

    opened_file_with_sim = open(sim_pickleTraffic, 'rb')
    data_from_pickle_sim = pickle.load(opened_file_with_sim)

    for range in data_from_pickle_sim.keys():
        plot_real_data(data_from_pickle_sim[range],
                   "/Users/natalia/Desktop/plots/heatmaps/new/moverate"+current_moverate+"/sim/" +range+"_heatmap.png",
                   title=range)


    #real_time_data = get_real_time_data("/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/fullBatch9",
    #                  "5", "Hunters Tryst", "9", 303, ['36234384'],['36234347'],['36237837'],['36247875'])


    #print(real_time_data[0][datetime.date(2018, 3, 8)]['6826'])

    real_time_data_pickle_path = "/Users/natalia/pickled_journey_times_realdata_9new.pickle"

    #process_real_time_data(real_time_data, all_pairs, real_time_data_pickle_path)
    average_journey_times_real_data = read_pickled_data_realtime(real_time_data_pickle_path)

    print(average_journey_times_real_data.keys())

    #range_str = '870.0-930.0'

    for range_str in average_journey_times_real_data.keys():
        #plot_real_data(average_journey_times_real_data[range_str], "/Users/natalia/Desktop/plots/heatmaps/newdata9_newnew/"+ range_str + "_heatmap.png", title=range_str)
        print("Finished " + range_str)

    difference_per_range  = dict()

    for range_str in average_journey_times_real_data.keys():
        reald = average_journey_times_real_data[range_str]
        simd = data_from_pickle_sim[range_str]

        difference_per_range[range_str] = get_diff_sum(reald, simd)
        diff = difference_between_heatmaps(reald, simd)

        plot_diff(diff, reald, simd,
                       "/Users/natalia/Desktop/plots/heatmaps/new/moverate"+current_moverate+"/diff/" + range_str + "_heatmap.png",
                       title=range_str)

    #plot_diff_per_range(difference_per_range, "/Users/natalia/Desktop/plots/heatmaps/new/moverate"+current_moverate+"/diffs.png", "/Users/natalia/Desktop/plots/heatmaps/new/moverate"+current_moverate+"/diffs.txt")

    #testPlotting()
    #data[busStart][busStop] = journeyTime
    #if busStart = busStop then journeyTime = 0
    #if busStart > busStop then journeyTime = 0






if __name__ == "__main__":
    main()
