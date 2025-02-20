export AIRFLOW_HOME=~/airflow;
mkdir -p $AIRFLOW_HOME/dags;
cp ./shell-tracer-autoinstrumentation/dags/rnaseq_dag.py $AIRFLOW_HOME/dags/;
airflow users create --username tracer --firstname Admin --lastname User --role Admin --email hello@tracer.bio;
airflow db migrate;
export AIRFLOW__CORE__LOAD_EXAMPLES=False;
airflow scheduler &;
airflow webserver --port 8080