import requests
import json
from itertools import groupby
import grein_loader as loader  
import time

# test grein request
#https://reactome.org/GSAServer/0.1/data/load/grein

list_grein_datasets = loader.load_overview(50)  
print(list_grein_datasets)

LOG_FILE_SUCCESS = "WHITELIST_RESULT_LOG.txt"  # successfull analysis 
LOGGING_SUCCESS = open(LOG_FILE_SUCCESS, 'w')

COUNTER = 0

# first testing stage request/load dataset

for grein_dataset in list_grein_datasets:
    
    print("COUNTER: " + str(COUNTER))

    if COUNTER > 5:
        print("Sleeping for 5 seconds")
        time.sleep(5)
        COUNTER = 0

    print("Testing dataset: " + grein_dataset["geo_accession"])
    url = 'https://reactome.org/GSAServer/0.1/data/load/grein'
    payload = '[{"name":"dataset_id","value":"' + grein_dataset["geo_accession"] + '"}]'


    headers = {
        'Content-Type': 'application/json'
    }
    time.sleep(4)
    response = requests.post(url, data=payload, headers=headers)

    # Check the response status code
    if response.status_code == 200:
        # Request was successful
        print('Response successful, :', response.text)
    else:
        # Request encountered an error
        print('Response unsuccessful:', response.text)
        COUNTER += 1
        continue
    

    url = 'https://reactome.org/GSAServer/0.1/data/summary/' + grein_dataset["geo_accession"]
    payload = '[{"name":"dataset_id","value":"' + grein_dataset["geo_accession"] + '"}]'
    time.sleep(4)
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        # Request was successful
        print('Summary Response successful:', response.text)
    else:
        # Request encountered an error
        print('Summary Response unsuccessful:', response.text)
        COUNTER += 1
        continue

    url = 'https://reactome.org/GSAServer/0.1/analysis'
    
    # check if samples and metadata list is the same size: 
    response_json = json.loads(response.text)
    metadata = response_json["sample_metadata"]
    list_valuse = metadata[0]["values"]
    
    
    if len(list_valuse) < 2: 
        print("Not enough metadata for analysis")
        COUNTER += 1
        continue

    groups = metadata[0]["values"]

    groups.sort()
    grouped_data = {key: list(group) for key, group in groupby(groups)}
    result = list(grouped_data.values()) # result is list of lists 

    if len(result) > 1:
        if len(result[0]) != len(result[1]):   # this appproach is not suitable for all datasets e.g. metadata list only contains same variables
            print("not the same size of groups")
            COUNTER += 1
            continue

    if len(result) < 1: 
        print("Not enough metadata for analysis")
        continue

    try:
        group1 = result[0][1]
        group2 = result[1][0]
        print("group1 exists:", group1)
        print("group2 exists:", group2)
    except IndexError:
        print("Not enough metadata for analysis")
        continue

    # Testing with different algorithms
    algorithms = ["PADOG","SSGSEA","CAMERA"]
    for algorithm in algorithms:
        payload = '{"methodName":"'+algorithm+'","datasets":[{"name":"test","type":"'+str(response_json["type"])+'","data":"'+str(response_json["id"])+'","design":{"analysisGroup":'+ str(metadata[0]["values"])+',"samples":'+ str(response_json["sample_ids"])+',"comparison":{"group1":"'+ str(group1) +'","group2":"'+ str(group2)+'"}}}],"parameters":[{"name":"use_interactors","value":"false"},{"name":"include_disease_pathways","value":"true"},{"name":"max_missing_values","value":"0.5"},{"name":"sample_groups","value":""},{"name":"discrete_norm_function","value":"TMM"},{"name":"continuous_norm_function","value":"none"},{"name":"create_reactome_visualization","value":"true"},{"name":"create_reports","value":"true"},{"name":"email","value":""},{"name":"reactome_server","value":"production"}]}'
        payload_new = payload.replace('\\', '')
        payload_new = payload_new.replace(' ', '')
        payload_new = payload_new.replace('\'', '\"')


        headers = {'Content-Type': 'application/json'}
        time.sleep(4)
        response = requests.post(url, data=payload_new, headers=headers)
        
        COUNTER += 1
        if response.status_code == 200:
            # Request was successful
            print('Response:, Analysis using '+algorithm, response.text)
            LOGGING_SUCCESS.write(grein_dataset["geo_accession"] + ";" + algorithm + "\n")
        else:
            # Request encountered an error
            print('Response Analysis unsuccessful using '+algorithm, response.text)
            continue
        
