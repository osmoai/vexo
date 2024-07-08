#!/bin/bash
set -xo pipefail

script_dir=$(dirname "$(realpath "$0")")
cd "$script_dir" || exit

cleanup() {
    echo "Cleaning up..."
}

# Trap the EXIT signal to call the cleanup function on script exit
trap cleanup EXIT

# Check if the dataset exists
dataset_exists=$(bq ls -d | grep -w "vexo")

# Create the dataset if it does not exist
if [ -z "$dataset_exists" ]; then
    bq mk "vexo"
    echo "Dataset 'vexo' created."
else
    echo "Dataset 'vexo' already exists."
fi

## Create connection

## Check if connection already exists

bq show --location=US --format=prettyjson --connection "vexo-connection" > /dev/null 2>&1

status=$?

## if the connection exists, continue. otherwise, create the connection

if [ $status -eq 0 ]
then
    echo "Connection vexo-connect already exists"
else
    echo "Creating connection vexo-connect"
    bq mk --connection --display_name="Vexo Connection" --connection_type=CLOUD_RESOURCE --location=US "vexo-connection"
fi

## Get Service Account associated to Connection

SERVICE_ACCOUNT=$(bq show --location=US --format=prettyjson --connection "vexo-connection" | jq -r '.cloudResource.serviceAccountId')

echo "Connection vexo-connect service account: ${SERVICE_ACCOUNT}"

## Give service account the cloud run invoker role (necessary for cloud functions gen2)

PROJ=$(gcloud config list --format 'value(core.project)')

####### Begin creating cloud functions

PERM="roles/cloudfunctions.invoker"

TIMEOUT=600s
MEMORY=512MB
MAX_INSTANCES=1000

while IFS=$'\t' read -r fn_call
do
    fn_name="vexo-${fn_call//_/-}"
    gcloud beta functions deploy $fn_name \
    --gen2 --region "us-east1" --entry-point $fn_call --runtime python311 --allow-unauthenticated --trigger-http \
    --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
    --update-labels package=vexo --update-labels function_type=remote_function --update-labels software_package=main
done < server_functions.csv

while IFS=$'\t' read -r fn_call in_args out_args
do
    fn_name="vexo-${fn_call//_/-}"
    echo "Creating function \"$fn_name\" with call \"$fn_call\" and input args \"$in_args\" and output args \"$out_args\""

    gcloud beta functions deploy $fn_name \
    --quiet --gen2 --region "us-east1" --entry-point $fn_call --runtime python311 --trigger-http \
    --memory=$MEMORY --timeout=$TIMEOUT --max-instances=$MAX_INSTANCES  \
    --update-labels package=vexo --update-labels function_type=remote_function --update-labels software_package=main

    CLOUD_TRIGGER_URL=$(gcloud beta functions describe $fn_name --gen2 --region "us-east1" --format=json | jq -r '.serviceConfig.uri')

    gcloud beta functions add-iam-policy-binding "$fn_name" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role=${PERM} --gen2

    gcloud run services add-iam-policy-binding "$fn_name" --region "us-east1" --member=serviceAccount:${SERVICE_ACCOUNT} --role="roles/run.invoker"

    bq query --use_legacy_sql=false --parameter="url::${CLOUD_TRIGGER_URL}" 'CREATE or REPLACE FUNCTION vexo.'$fn_call'('$in_args') RETURNS '$out_args' REMOTE WITH CONNECTION `us.vexo-connection` OPTIONS (endpoint = @url, max_batching_rows = 2500)'
done < bigquery_functions.csv