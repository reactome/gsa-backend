import click
import os
import requests
import grein_loader.exceptions
import grein_loader as loader
import logging

LOGGER = logging.getLogger(__name__)


@click.command()
@click.option("--path_blacklist", default=None, help="Path to the blacklist file")
@click.option("--path_whitelist", default=None, help="Path to the whitelist file")
def update_lists(path_blacklist, path_whitelist):
    if not path_blacklist:
        path_blacklist = os.getenv("SEARCH_INDEX_BLACKLIST", "../")

    if not path_whitelist:
        path_whitelist = os.getenv("SEARCH_INDEX_WHITELIST", "../")

    grein_datasets = loader.load_overview()
    blacklist = open(path_blacklist, "r+")
    whitelist = open(path_whitelist, "r+")

    grein_datasets_ = [d["geo_accession"] for d in grein_datasets]

    checked_datasets = set(line.strip() for line in whitelist)
    checked_datasets.update(line.strip() for line in blacklist)
    datasets_to_check = list(set(grein_datasets_) - checked_datasets)

    counter_blacklist = 0
    counter_whitelist = 0
    if len(datasets_to_check) == 0:
        LOGGER.info("No new datasets from GREIN")
    else:
        print(str(len(datasets_to_check)) + " new datasets")
        LOGGER.info(str(len(datasets_to_check)) + " new datasets")
        for dataset_id in datasets_to_check:
            try:
                des, metadata, cm = loader.load_dataset(dataset_id)
                if metadata == "":
                    blacklist.write(dataset_id + "\n")
                    counter_blacklist += 1
                else:
                    no_metadata_keys = len(metadata.keys())
                    no_cm_keys = len(cm.axes[1]) - 2
                    if des["Summary"] == "":
                        blacklist.write(dataset_id + "\n")
                        counter_blacklist += 1
                    else:
                        if no_cm_keys == no_metadata_keys:
                            whitelist.write(dataset_id + "\n")
                            counter_whitelist += 1
                        else:
                            blacklist.write(dataset_id + "\n")
                            counter_blacklist += 1
            except requests.exceptions.HTTPError as e:
                LOGGER.error("Error fetching data from GREIN", e)
                blacklist.write(dataset_id + "\n")
            except grein_loader.exceptions.GreinLoaderException as eg:
                LOGGER.error("Error fetching data from GREIN via grein_loader", eg)
                blacklist.write(dataset_id+"\n")
    LOGGER.info("Lists updated; Whitelist: " + str(counter_whitelist) + " Blacklist: " + str(counter_blacklist))
    blacklist.close()
    whitelist.close()
