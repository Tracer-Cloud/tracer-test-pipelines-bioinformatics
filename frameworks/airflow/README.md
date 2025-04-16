# Airflow Testing For Tracer 
- We assume tracer runs in the background
- When testing on EC2 AWS you need to be logged in as an ubuntu user and execute all commands with sudo 


## Instructions To Run Airflow With Docker
- You need to use docker

```bash
cd airflow &&
sudo docker-compose up airflow-init &&
sudo docker compose up -d    
```

View docker process 
```bash
cd airflow-docker &&
 sudo docker compose ps
```

```bash
# Restart up and down to apply changes
sudo docker compose down
sudo docker compose up -d
```

## Run the pipeline 
```bash
tracer init --pipeline-name airflow_test --environment sandbox_test --user-operator vincent --pipeline-type rnaseq
```

## How to find the results of a specific pipeline run

1. Unpause the pipeline that you want to run
```bash
sudo docker compose run --rm airflow-cli dags unpause pipeline_rnaseq_2
```

2. Run that particular pipeline
```bash
sudo docker compose run --rm airflow-cli dags trigger pipeline_rnaseq_2 --run-id=my_custom_run_001
```
3. Get the results of that particular pipeline run
```bash
sudo docker compose run --rm airflow-cli tasks states-for-dag-run pipeline_rnaseq_2 my_custom_run_001
```

## Running Airflow Tasks
```bash
# Reapply changes if you need to trobleshoot
sudo docker compose restart airflow-scheduler
```

```bash
# Unpause the DAG
sudo docker compose run --rm airflow-cli dags airflow dags unpause pipeline_rnaseq_2

# Trigger the DAG (will reply false)
sudo docker compose run --rm airflow-cli dags trigger pipeline_rnaseq_2

# Check the status
sudo docker compose run --rm airflow-cli dags state pipeline_rnaseq_2

# List runs
sudo docker compose run --rm airflow-cli dags list-runs --dag-id pipeline_rnaseq_2
```

## Clean up old pipeline runs
```bash
sudo docker compose run --rm airflow-cli airflow dags delete pipeline_rnaseq_2
```


## Clean up remove logs

```bash
rm -rf logs/* 
```

## Webserver is exposed at:
http://127.0.0.1:8080/login/

Password and username:
- Username: airflow
- Password: airflow


# Requirements and installations:
brew install jq
