#!/bin/bash

#-----------------------------------------------------
# This script is used to build the Docker images used
# for the ReactomeGSA Analysis System backend
#
# The script is interactive with user input expected
# at different levels.
#
# The script is tested under Ubuntu 19.04
# ----------------------------------------------------

# CONFIG

# The Docker Hub username to tag the images with
DOCKER_USER="jgriss"

# make sure the working directory is the "docker" directory
if [ ! -d "../reactome_analysis_api" -o ! -d "../reactome_analysis_utils" -o ! -d "../reactome_analysis_worker" ]; then
    echo "Error: Failed to find required python directories. The script must be executed from within the 'docker' directory."
    exit 1
fi

CUR_DIR=`pwd`

# general function to check for errors
function check_error {
	RET="$1"
	MSG="$2"

	if [ ${RET} != 0 ]; then
		echo "Error: ${MSG}"
		exit 1
	fi
}

# Get the complete user input
echo -n "Update Reactome mappings [y/N]: "
read UPDATE_REACTOME

echo -n "Rebuild python packages [Y/n]: "
read UPDATE_PYTHON

echo -n "Rebuild public API [version/N]: "
read REBUILD_API

echo -n "Rebuild worker [version/N]: "
read REBUILD_WORKER

echo -n "Rebuild report [version/N]: "
read REBUILD_REPORT

# get sudo privileges required to work with docker
sudo echo ""

# Update the Reactome mappings (if set)
# This function downloads all required files from the Reactome
# content service
if [ "${UPDATE_REACTOME}" == "y" ]; then
	echo "Updating REACTOME mappings..."

	wget -O "UniProt2Reactome_All_Levels.txt" "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"
	check_error "$?" "Failed to get UniProt mappings"

	wget -O "Ensembl2Reactome_All_Levels.txt" "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt"
	check_error "$?" "Failed to get ENSEMBL mappings"

	wget -O "IntAct_Static.txt" "https://reactome.org/download/current/IntAct_Static.txt"
	check_erro "$?" "Failed to get IntAct mappings"
fi

# Make sure all files are present
for MAPPING_FILE in "UniProt2Reactome_All_Levels.txt" "Ensembl2Reactome_All_Levels.txt" "IntAct_Static.txt"; do
    if [ ! -e "${MAPPING_FILE}" ]; then
        echo "Error: Mapping file '${MAPPING_FILE}' does not exist"
        exit 1
    fi
done

# rebuild the python packages if set
if [ "${UPDATE_PYTHON}" != "n" ]; then
	echo "Updating python packages..."

	for PROJECT in "reactome_analysis_api" "reactome_analysis_utils" "reactome_analysis_worker" "reactome_analysis_report"; do
		echo -n "  ${PROJECT}..."
		cd "../${PROJECT}"
        check_error $? "Failed to access directory"

        # remove previous builds
        if [ -d "dist" ]; then
            rm -rf "dist"
        fi

        # re-build the package
		python3 setup.py bdist_wheel &> /dev/null
		check_error $? "Failed"
		echo "OK"

        cd "${CUR_DIR}"
	done
fi

# copy the python packages to the docker build context
# For Docker, all files required within the images have
# to be in the same folder as the Docker file
find ../ -name "*.whl" -exec cp {} ${CUR_DIR}/ \;
check_error $? "Failed to copy .whl files"

cp ../reactome_analysis_worker/install_libraries.R ${CUR_DIR}/
check_error $? "Failed to copy install_libraries.R"

cp ../reactome_analysis_report/install_libraries.R ${CUR_DIR}/install_report_libraries.R
check_error $? "Failed to copy install_libraries.R for report"

cp ../reactome_analysis_report/texlive.profile ${CUR_DIR}/

# Rebuild the images
if [ -n "${REBUILD_API}" -a "${REBUILD_API}" != "n" ]; then
	sudo docker build --no-cache -t ${DOCKER_USER}/reactome-analysis_public-api:${REBUILD_API} -f Dockerfile.public_api .
	check_error $? "Failed to build public API"
fi

if [ -n "${REBUILD_WORKER}" -a "${REBUILD_WORKER}" != "n" ]; then
	sudo docker build -t ${DOCKER_USER}/reactome-analysis_worker:${REBUILD_WORKER} -f Dockerfile.worker .
	check_error $? "Failed to build worker"
fi

if [ -n "${REBUILD_REPORT}" -a "${REBUILD_REPORT}" != "n" ]; then
    # TODO: Fix worker version reference
    if [ -z "${REBUILD_WORKER}" -o "${REBUILD_WORKER}" == "n" ]; then
        echo -n "Enter reactome_analysis_worker image version: "
        read WORKER_VERSION
    else
        WORKER_VERSION="${REBUILD_WORKER}"
    fi

    # replace the worker version in the Docker file
    sed "s/THE_WORKER_VERSION/${WORKER_VERSION}/" Dockerfile.report > Dockerfile.report_version

	sudo docker build -t ${DOCKER_USER}/reactome-analysis_report:${REBUILD_REPORT} -f Dockerfile.report_version .

    check_error $? "Failed to build report"

    # remove the versioned file
    rm Dockerfile.report_version	
fi