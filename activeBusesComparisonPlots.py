__author__ = 'natalia'

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


# Number of active buses vs time of the day
# Plotted for four cases:
# simulation data(obtained through carma interface),        carmadata
# simulation data (obtained from java_hacking measures),    rawjavadata
# timetable data (Lotian Api),                              timetabledata
# real-time data (Lothian Api)                              realtimedata
#


# returns [times, values]
def get_carma_data(path_to_csv_file_from_carma):
    opened_data_file = open(path_to_csv_file_from_carma)
    lines_of_data = opened_data_file.readlines()
    # time;value;error1;error2;
    # ERROR1 is the standard deviation. ERROR2 is the standard error of the mean and equal to std / sqrt(n)
    # where n is the sample size.

    list_of_times = list()
    list_of_means = list()
    list_of_stds = list()
    list_of_sems = list()

    for lineofdata in lines_of_data:
        #print(lineofdata)
        splat = lineofdata.split(";")
        time_of_simulation = splat[0]
        active_buses_mean = splat[1]
        active_buses_std = splat[2]
        active_buses_sem = splat[3]

        list_of_times.append(float(time_of_simulation))
        list_of_means.append(float(active_buses_mean))
        list_of_stds.append(float(active_buses_std))
        list_of_sems.append(float(active_buses_sem))

    return {"simulation_times": list_of_times, "activeBusesMeans": list_of_means, "activeBusesStds": list_of_stds,
            "activeBusesSems": list_of_sems}


def test_carmadata():
    pathto_carmadata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/" \
                       "carmadata/busActive16_02_2018__120551.csv"
    get_carma_data(pathto_carmadata)


# returns [times, values]
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


def test_rawjavadata():
    pathtoRawjavadata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/rawjavadata/2018.02.16.12.05_moveUpdates.txt"
    get_rawjavadata(pathtoRawjavadata)


# returns [times, values]
def get_timetable_data(path_to_timetable_data):
    print(path_to_timetable_data)
    return []


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


# creates a pickled and reorganized version of live vehicle location data obtained from lothian api
# outputs one pickle file per service
# run only once for each batch of data, then use the pickled versions
def reorganize_live_data_collected(path_to_data, dataset_name):
    print(path_to_data)
    month_number = dict((v, k) for k, v in enumerate(calendar.month_abbr))
    files = os.listdir(path_to_data)
    valid_filename_regexp = re.compile('dataset9_[0-9]{4}-[A-z]{3}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}\.txt')

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
                "fullBatch9/"+dataset_name+"_service" + "_" + serviceName + ".pickle",
                "wb")
            pickle.dump(reorganized[serviceName], pickle_out)
            pickle_out.close()
        except TypeError:
            print("Type error for service ")
            continue

    return


# returns [times, values]
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
                    #if(last_seen_snapshot[1]['next_stop_id'] in last_stop_ids):
                        #data is valid -> carry on processing

                        first_seen = first_seen_snapshot[0]
                        last_seen = last_seen_snapshot[0]

                        #For fist seen, if next stop id is the initial stop, choose the last timesnap at with that value of next_stop_id
                        if first_seen_snapshot[1]['next_stop_id'] in first_stop_ids:
                            snapshots_with_first_stop_id = list()
                            for snapshot in data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id]:
                                if snapshot[1]['next_stop_id'] in first_stop_ids:
                                    snapshots_with_first_stop_id.append(snapshot)
                                    first_seen = max(snapshots_with_first_stop_id,
                        key=operator.itemgetter(0))[0]


                        #For fist seen, if next stop id is the second stop, choose the earliest timesnap at with that value of next_stop_id
                        elif first_seen_snapshot[1]['next_stop_id'] in second_stop_ids:
                            first_seen = first_seen_snapshot[0]

                        else:
                            raise Exception("Illegal State Exception")

                        #For last seen, choose the earliest timesnap at with the value of next_stop_id == last_stop_id
                        # snaps_with_last_stop_id = list()
                        # for snapshot in data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id]:
                        #     if snapshot[1]['next_stop_id'] in last_stop_ids:
                        #         snaps_with_last_stop_id.append(snapshot)
                        # last_seen = min(snaps_with_last_stop_id, key=operator.itemgetter(0))[0]
                        #

                        #last_seen = last_seen_snapshot[0]

                        #first_seen = first_seen_snapshot[0]

                        #print("First seen " + str(first_seen) + ", last seen " + str(last_seen))

                        first_and_last_seen[type_of_day][date_value][journey_id]['first_seen'] = first_seen
                        first_and_last_seen[type_of_day][date_value][journey_id]['last_seen'] = last_seen
                    # else:
                    #     first_and_last_seen[type_of_day][date_value][journey_id]['first_seen'] = None
                    #     first_and_last_seen[type_of_day][date_value][journey_id]['last_seen'] = None
                    #     continue
                elif (first_seen_snapshot[1]['next_stop_id'] in wrong_stop_ids or first_seen_snapshot[1]['next_stop_id']in wrong_stop_ids or last_seen_snapshot[1]['next_stop_id'] in wrong_stop_ids):
                    print("Nonexistant stop bus " + str(journey_id) + " first " + str(
                       first_seen_snapshot[1]['next_stop_id']) + " last " + last_seen_snapshot[1]['next_stop_id'])
                    first_seen = None
                    last_seen = None
                    first_and_last_seen[type_of_day][date_value][journey_id]['first_seen'] = first_seen
                    first_and_last_seen[type_of_day][date_value][journey_id]['last_seen'] = last_seen
                    continue
                else:
                    print("Misbehaved bus " + str(journey_id) + " first " + str(first_seen_snapshot[1]['next_stop_id']) + " last " + last_seen_snapshot[1]['next_stop_id'])
                    first_seen = first_seen_snapshot[0]
                    last_seen = last_seen_snapshot[0]
                    #first_seen = None
                    #last_seen = None
                    first_and_last_seen[type_of_day][date_value][journey_id]['first_seen'] = first_seen
                    first_and_last_seen[type_of_day][date_value][journey_id]['last_seen'] = last_seen
                    continue


    


    list_of_simulation_timestep_values = list()
    for i in range(0, 1440, 1):
        list_of_simulation_timestep_values.append(i)

    active_counts_per_daytype_per_day_per_time = defaultdict(d1)
    inactive_counts_per_daytype_per_day_per_time = defaultdict(d1)

    journeys_active_at_night = set()

    for type_of_day in sorted(data_by_daytype_by_day_by_journey.keys()):
        for date_value in sorted(data_by_daytype_by_day_by_journey[type_of_day].keys()):
            active_buses_per_min = dict()
            inactive_buses_per_min = dict()
            for i in range(0, 1440, 1):
                active_count = 0
                inactive_count = 0
                for journey_id in sorted(first_and_last_seen[type_of_day][date_value]):
                    #print("!!! " + str(data_by_daytype_by_day_by_journey[type_of_day][date_value][journey_id]))
                    first_seen =  first_and_last_seen[type_of_day][date_value][journey_id]['first_seen']
                    last_seen =   first_and_last_seen[type_of_day][date_value][journey_id]['last_seen']

                    #print("F seen" + str(first_seen))
                    #print("L seen" + str(last_seen))

                    if(first_seen == None or last_seen == None):
                        continue

                    first_seen_min = first_seen.hour * 60 + first_seen.minute
                    last_seen_min = last_seen.hour * 60 + last_seen.minute
                    if (first_seen_min < i < last_seen_min):

                        if (0.0 < i < timetable_start_in_min):
                            #print("Journey " + str(journey_id) + " active at " + str(i/60.0) )
                            journeys_active_at_night.add(journey_id)

                        else:
                            active_count += 1

                    else:
                        inactive_count += 1


                active_buses_per_min[i] = active_count
                inactive_buses_per_min[i] = inactive_count
            active_counts_per_daytype_per_day_per_time[type_of_day][date_value] = active_buses_per_min
            inactive_counts_per_daytype_per_day_per_time[type_of_day][date_value] = inactive_buses_per_min
    print("Journeys active out of timetable: " + str(journeys_active_at_night))

    active_buses_per_day_type = dict()
    inactive_buses_per_day_type = dict()

    for day_type in active_counts_per_daytype_per_day_per_time.keys():
        lists_of_active_bus_counts_per_sim_time = defaultdict(list)
        lists_of_inactive_bus_counts_per_sim_time = defaultdict(list)
        for date_value in active_counts_per_daytype_per_day_per_time[day_type].keys():
            for min_val in active_counts_per_daytype_per_day_per_time[day_type][date_value].keys():
                lists_of_active_bus_counts_per_sim_time[min_val].append(
                    active_counts_per_daytype_per_day_per_time[day_type][date_value][min_val])
                lists_of_inactive_bus_counts_per_sim_time[min_val].append(
                    inactive_counts_per_daytype_per_day_per_time[day_type][date_value][min_val])
        active_buses_per_day_type[day_type] = lists_of_active_bus_counts_per_sim_time
        inactive_buses_per_day_type[day_type] = lists_of_inactive_bus_counts_per_sim_time

    list_of_active_vals_with_errors_per_daytype = defaultdict(list)
    list_of_active_means_per_daytype = defaultdict(list)
    list_of_active_stds_per_daytype = defaultdict(list)
    list_of_active_sems_per_daytype = defaultdict(list)

    list_of_inactive_vals_with_errors_per_daytype = defaultdict(list)
    list_of_inactive_means_per_daytype = defaultdict(list)
    list_of_inactive_stds_per_daytype = defaultdict(list)
    list_of_inactive_sems_per_daytype = defaultdict(list)

    for day_type in active_buses_per_day_type.keys():
        for simul_time in sorted(active_buses_per_day_type[day_type].keys()):
            mean = np.mean(active_buses_per_day_type[day_type][simul_time])
            #print("? " + str(day_type) + " " + str(simul_time/60.0) + " " + str(active_buses_per_day_type[day_type][simul_time]))
            std = np.std(active_buses_per_day_type[day_type][simul_time], ddof=0)
            semval = stats.sem(active_buses_per_day_type[day_type][simul_time], ddof=0)
            list_of_active_vals_with_errors_per_daytype[day_type].append([mean, std, semval])
            list_of_active_means_per_daytype[day_type].append(mean)
            list_of_active_stds_per_daytype[day_type].append(std)
            list_of_active_sems_per_daytype[day_type].append(semval)

    for day_type in inactive_buses_per_day_type.keys():
        for simul_time in sorted(inactive_buses_per_day_type[day_type].keys()):
            mean = np.mean(inactive_buses_per_day_type[day_type][simul_time])

            std = np.std(inactive_buses_per_day_type[day_type][simul_time], ddof=0)
            semval = stats.sem(inactive_buses_per_day_type[day_type][simul_time], ddof=0)
            list_of_inactive_vals_with_errors_per_daytype[day_type].append([mean, std, semval])
            list_of_inactive_means_per_daytype[day_type].append(mean)
            list_of_inactive_stds_per_daytype[day_type].append(std)
            list_of_inactive_sems_per_daytype[day_type].append(semval)

    print("Returning from get_real_time_data ... There are " + str(len(active_counts_per_daytype_per_day_per_time[0])) + " weekday days of data ")



    return {"simulation_times": list_of_simulation_timestep_values,
            "activeBusesMeans": list_of_active_means_per_daytype[0],
            "activeBusesStds": list_of_active_stds_per_daytype[0],
            "activeBusesSems": list_of_active_sems_per_daytype[0],
            "inactiveBusesMeans": list_of_inactive_means_per_daytype[0],
            "inactiveBusesStds": list_of_inactive_stds_per_daytype[0],
            "inactiveBusesSems": list_of_inactive_sems_per_daytype[0]}


def test_get_real_time_data():
    # reorganize_live_data_collected("/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/set4")
    get_real_time_data("/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/pickled_new_version",
                       "5", "Hunters Tryst")


def get_data_from_sources():
    return
    # pathto_carmadata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/carmadata/" \
    #                   "busActive16_02_2018__120551.csv"
    # carmadata_loaded = get_carma_data(pathto_carmadata)
    # pathto_rawjavadata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/rawjavadata/" \
    #                    "2018.02.16.12.05_moveUpdates.txt"
    # raw_java_data_loaded = get_rawjavadata(pathto_rawjavadata)
    # pathto_timetabledata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/" \
    #                       "timetabledata/timetables.p"  # pickled
    # timetable_data_loaded = get_timetable_data(pathto_timetabledata)
    # pathto_realtimedata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/set4"  # folder
    # realtimedata_loaded = get_real_time_data(pathto_realtimedata)
    # pathto_pickled_realtimedata = "/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/" \
    #                              "realtimedata/reorganized_pickled"  # folder, one pickle file per service


def rgb_frac(rgb_colour_255):
    red = float(rgb_colour_255[0]) / 255.0
    green = float(rgb_colour_255[1]) / 255.0
    blue = float(rgb_colour_255[2]) / 255.0
    opac = float(rgb_colour_255[3])
    return (red, green, blue, opac)

from matplotlib import colors as mcolors

def get_colour(name, transparency):
    mycolor=mcolors.to_rgba(name)
    new_color = rgb_frac((mycolor[0]*255, mycolor[1]*255, mycolor[2]*255, transparency))
    return new_color

from scipy.interpolate import spline
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib
def plot_active_and_inactive(datasets, names):



    font = {'family': 'normal',
            'sans-serif':'Helvetica',
            'weight': 'medium',
            'size': 16}

    matplotlib.rc('font', **font)

    colours_per_index = dict()
    colours_per_index[0] = {
        'solid': [get_colour('b', 1.0), get_colour('blueviolet', 1.0), get_colour('midnightblue', 1.0)],
        'transparent': [get_colour('b', 0.15), get_colour('blueviolet', 0.2), get_colour('midnightblue', 0.08)]}
    colours_per_index[1] = {
        'solid': [get_colour('red', 1.0), get_colour('coral', 1.0), get_colour('firebrick', 1.0)],
        'transparent': [get_colour('red', 0.1), get_colour('coral', 0.2), get_colour('firebrick', 0.08)]}
    colours_per_index[2] = {
        'solid': [get_colour('darkgreen', 1.0), get_colour('yellow', 1.0), get_colour('darkgreen', 1.0)],
        'transparent': [get_colour('darkgreen', 0.2), get_colour('yellow', 0.2), get_colour('darkgreen', 0.08)]}

    line_kinds =[':', '--', '-']
    patches  = []
    lines= []
    for index in range(len(datasets)):
        print(names[index])
        print(datasets[index])

        # for value in datasets[index]["simulation_times"]:
        #     openedFileSimTimes.write(str(value) + "\n")
        # for value in datasets[index]["activeBusesMeans"]:
        #     openedFileActiveSIm.write(str(value) + "\n")
        # for value in datasets[index]["activeBusesStds"]:
        #     openedFileActiveStd.write(str(value) + "\n")
        # for value in datasets[index]["activeBusesSems"]:
        #     openedFileActiveSem.write(str(value) + "\n")

        x_axis_vals = np.array(datasets[index]["simulation_times"]) / 60.0
        y_vals_active_buses_means = np.array(datasets[index]["activeBusesMeans"])
        y_vals_active_buses_stds = np.array(datasets[index]["activeBusesStds"])
        y_vals_active_buses_sems = np.array(datasets[index]["activeBusesSems"])

        # x_smooth = np.linspace(x_axis_vals.min(), x_axis_vals.max(), 96)
        # y_smooth1= spline(x_axis_vals, y_vals_active_buses_means - y_vals_active_buses_stds, x_smooth)
        # y_smooth2= spline(x_axis_vals,y_vals_active_buses_stds + y_vals_active_buses_means , x_smooth)

        plt.plot(x_axis_vals, y_vals_active_buses_means - y_vals_active_buses_stds, line_kinds[index], label=names[index] + " std",
                 linewidth=1.0, color=colours_per_index[index]['transparent'][0])

        plt.plot(x_axis_vals, y_vals_active_buses_means + y_vals_active_buses_stds,line_kinds[index],
                 label=names[index] + " std", linewidth=1.0, color=colours_per_index[index]['transparent'][0])

        # plt.plot(x_smooth, y_smooth1, line_kinds[index], label=names[index] + " std",
        #          linewidth=1.0, color=colours_per_index[index]['transparent'][0])
        #
        #
        # plt.plot(x_smooth, y_smooth2, line_kinds[index], label=names[index] + " std",
        #          linewidth=1.0, color=colours_per_index[index]['transparent'][0])

        # plt.fill_between(x_smooth, y_smooth1,
        #                  y_smooth2,
        #                  color=colours_per_index[index]['transparent'][0])

        plt.fill_between(x_axis_vals, y_vals_active_buses_means - y_vals_active_buses_stds,
                         y_vals_active_buses_means + y_vals_active_buses_stds,
                         color=colours_per_index[index]['transparent'][0])

        legend_patch = mpatches.Patch(color=colours_per_index[index]['transparent'][0], label=names[index] + " std")
        patches.append(legend_patch)

        # plt.plot(x_axis_vals, y_vals_active_buses_means - y_vals_active_buses_sems, '--', label=names[index] + " sem",
        #          linewidth=0.5, color=colours_per_index[index]['solid'][1])
        # plt.plot(x_axis_vals, y_vals_active_buses_means + y_vals_active_buses_sems, '--', label=names[index] + " sem",
        #          linewidth=0.5, color=colours_per_index[index]['solid'][1])
        # plt.fill_between(x_axis_vals, y_vals_active_buses_means - y_vals_active_buses_sems,
        #                  y_vals_active_buses_means + y_vals_active_buses_sems,
        #                  color=colours_per_index[index]['transparent'][1])

        # plt.plot(datasets[index]["simulation_times"], datasets[index]["activeBusesStds"], label="Active buses means" )
        # plt.plot(datasets[index]["simulation_times"], datasets[index]["activeBusesStds"], label="Active buses means" )

        # plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')

    for index in range(len(datasets)):
        x_axis_vals = np.array(datasets[index]["simulation_times"]) / 60.0
        y_vals_active_buses_means = np.array(datasets[index]["activeBusesMeans"])

        #x_smooth = np.linspace(x_axis_vals.min(), x_axis_vals.max(), 96)
        #y_smooth = spline(x_axis_vals, y_vals_active_buses_means, x_smooth)

        y_vals_active_buses_stds = np.array(datasets[index]["activeBusesStds"])
        y_vals_active_buses_sems = np.array(datasets[index]["activeBusesSems"])
        plt.plot(x_axis_vals, y_vals_active_buses_means,line_kinds[index], label=names[index] + "",
                  color=colours_per_index[index]['solid'][2], linewidth=1.5)

        #plt.plot(x_smooth, y_smooth,line_kinds[index], label=names[index] + "",
        #                 color=colours_per_index[index]['solid'][2], linewidth=1.5)

        line = mlines.Line2D([],[], linestyle=line_kinds[index], color=colours_per_index[index]['solid'][2],linewidth=1.5, label=names[index] + "" )
        lines.append(line)



    plt.grid()




    plt.xlabel("Time of day [hour]")
    plt.ylabel("Number of active buses")

    lines_and_patches = []
    for idx in range(len(patches)):
        lines_and_patches.append(lines[idx])
        lines_and_patches.append(patches[idx])

    plt.legend(handles=lines_and_patches, fontsize="x-small", ncol=2)
    plt.axis([4.5, 15.5, 0, 10])
    plt.xticks(np.arange(5, 16, 1))


    plt.yticks(np.arange(0,10,1))

    #
    #plt.setp(ax3.yaxis.get_label(), visible=True, rotation=30, ha='right')

    plt.show()


def main():
    print("Running BusActive comparison plots script")
    # test_rawjavadata()
    pathto_rawjavadata = "/Users/natalia/Desktop/plots/results_files/brunstane_weighted/" \
                         "2018.03.08.18.22_moveUpdateD.txt"
    #
    raw_java_data = get_rawjavadata(pathto_rawjavadata)

    csv_data_with_traffic = get_carma_data("/Users/natalia/Work/version27Aug17/workspace/LothianBusesCarma/output/busActive26_03_2018__220744.csv")
    csv_data_without_traffic = get_carma_data("/Users/natalia/Work/version27Aug17/workspace/LothianBusesCarma/output/busActive26_03_2018__232525.csv")
    #reorganize_live_data_collected("/Users/natalia/Work/LothianBusesDataPaper/data/downloaded_and_processed/fromDell/liveDataGathered_newTimetable_Batch9/", "9")


    try:
        pickle_out = open(
            "/Users/natalia/Work/LothianBusesDataPaper/" + "temporary_rawjavadata2"+ ".pickle",
            "wb")
        pickle.dump(raw_java_data, pickle_out)
        pickle_out.close()
    except TypeError:
        print("Error pickling rawjavadata ")



    #

    #cut the data off at 5:03 (303 min)
    real_time_data = get_real_time_data("/Users/natalia/Work/LothianBusesDataPaper/busActivePlots/data/realtimedata/fullBatch9",
                      "5", "Hunters Tryst", "9", 303, ['36234384'],['36234347'],['36237837'],['36247875'])


    #
    #
    # try:
    #     pickle_out = open(
    #         "/Users/natalia/Work/LothianBusesDataPaper/" + "temporary_realtimedata"+ ".pickle",
    #         "wb")
    #     pickle.dump(real_time_data, pickle_out)
    #     pickle_out.close()
    # except TypeError:
    #   print("Error pickling realtimedata ")

    # opened_file_rawjavadata = open("/Users/natalia/Work/LothianBusesDataPaper/" + "temporary_rawjavadata" + ".pickle",
    #                                'rb')
    # loaded_rawjavadata_pickle = pickle.load(opened_file_rawjavadata)

    #opened_file_realtimedata = open("/Users/natalia/Work/LothianBusesDataPaper/" + "temporary_realtimedata" + ".pickle",
    #                                'rb')
    #loaded_realtime_pickle = pickle.load(opened_file_realtimedata)

    #main plotting function
    plot_active_and_inactive([csv_data_without_traffic,csv_data_with_traffic, real_time_data], ["Simul","Simul with traffic", "Data"])

    return


if __name__ == "__main__":
    main()
