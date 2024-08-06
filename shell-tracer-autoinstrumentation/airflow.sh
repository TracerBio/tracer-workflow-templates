#!/bin/bash

# Directory where Airflow data and configuration will be stored
PATH_DIR=/workspace/tracer-workflow-templates

# Directory where Airflow data and configuration will be stored
AIRFLOW_DIR=$PATH_DIR/airflow-docker

# Create the directory if it doesn't exist
mkdir -p $AIRFLOW_DIR

# Navigate to the Airflow directory
cd $AIRFLOW_DIR

# Download the official Apache Airflow docker-compose.yaml
curl -LfO 'https://airflow.apache.org/docs/apache-airflow/2.0.2/docker-compose.yaml'

# Create a .env file for environment variables
echo -e "AIRFLOW_UID=$(id -u)\nAIRFLOW_GID=0" > .env

# Initialize the Airflow database
docker-compose up airflow-init

# Start the Airflow services in detached mode
docker-compose up -d

# Provide feedback to the user
echo "Apache Airflow is set up and running."
echo "Access the Airflow web interface at http://localhost:8080"