#!/bin/bash
set -xo pipefail

script_dir=$(dirname "$(realpath "$0")")
cd "$script_dir" || exit

cleanup() {
  echo "Cleaning up..."
}

####### Begin creating cloud functions

PERM="roles/cloudfunctions.invoker"

TIMEOUT=600s
MEMORY=512MB
MAX_INSTANCES=1000

PROJ=$(gcloud config list --format 'value(core.project)')
   
gcloud beta functions deploy vexo-smiles-png \
    --gen2 --region "us-east1" --entry-point smiles_png --runtime python311 --allow-unauthenticated --trigger-http \
    --update-labels package=vexo --update-labels function_type=remote_function --update-labels software_package=main


gcloud beta functions deploy vexo-smart-substructure \
    --gen2 --region "us-east1" --entry-point smart_substructure --runtime python311 --allow-unauthenticated --trigger-http \
    --update-labels package=vexo --update-labels function_type=remote_function --update-labels software_package=main