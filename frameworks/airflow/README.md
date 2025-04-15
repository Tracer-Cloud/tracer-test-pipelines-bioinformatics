# Airflow Testing For Tracer 
- We assume tracer runs in the background


## How to find the results of a specific pipeline run
1. First run that particular pipeline
```bash
docker compose run --rm airflow-cli dags trigger pipeline_rnaseq --run-id=my_custom_run_001
```
2. Get the results of that particular pipeline run
```bash
docker compose run --rm airflow-cli tasks states-for-dag-run pipeline_rnaseq my_custom_run_001
```

## Instructions To Run Airflow With Docker
- You need to use docker

```bash
cd airflow-docker &&
docker-compose up airflow-init &&
docker compose up -d    
```

View docker process 
```bash
cd airflow-docker &&
 docker compose ps
```

```bash
# Restart up and down to apply changes
docker compose down
docker compose up -d
```


## Running Airflow Tasks
```bash
# Reapply changes if you need to trobleshoot
docker compose restart airflow-scheduler
```

```bash
# Unpause the DAG
docker compose run --rm airflow-cli dags airflow dags unpause pipeline_rnaseq

# Trigger the DAG (will reply false)
docker compose run --rm airflow-cli dags trigger pipeline_rnaseq

# Check the status
docker compose run --rm airflow-cli dags state pipeline_rnaseq

# List runs
docker compose run --rm airflow-cli dags list-runs --dag-id pipeline_rnaseq
```

## Clean up old pipeline runs
```bash
docker compose run --rm airflow-cli airflow dags delete pipeline_rnaseq
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
