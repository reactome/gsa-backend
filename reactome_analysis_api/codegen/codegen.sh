#!/bin/bash

####################################################
# This script is used to re-generate the code
# (mainly the data model) in the reactome_analysis_api
# based on the swagger definition
####################################################

# make sure the script is run from its own directory
if [ ! -f "../.swagger-codegen-ignore" ]; then
	echo "Error: The script must be executed from within its directory"
	exit 1
fi

# make sure java is installed
if [ ! -e "java" ]; then
	echo "Failed to find Java"
	exit 1
fi

# move to the parent directory
cd ..

CODEGEN_JAR="codegen/swagger-codegen-cli.jar"

if [ ! -f "${CODEGEN_JAR}" ]; then
	echo "Error: Swagger Codegen jar file not found."
	exit 1
fi

java -jar "${CODEGEN_JAR}" generate \
	-l python \
	--additional-properties "packageName=reactome_analysis_api, projectName=reactome_analysis-public_api, packageVersion=0.1, library=flask" \
	-i analysis_system-public_api.swagger.yaml \
	-l python-flask \
	--ignore-file-override "../.swagger-codegen-ignore"

