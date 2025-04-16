# Airflow Testing For Tracer 
- We assume tracer runs in the background
- When testing on EC2 AWS you need to be logged in as an ubuntu user and execute all commands with sudo 


## Instructions To Run Airflow With Docker
- You need to use docker

```bash
sudo docker compose up -d    &&
sudo docker-compose up airflow-init 
```

```bash
# View running docker process 
 sudo docker compose ps
```

```bash
# Ensure the database is running
sudo docker compose run --rm airflow-cli airflow db init
```

```bash
# Ensure the airflow scheduler is up
sudo docker compose up -d airflow-scheduler
```

```bash
# Restart up and down to apply changes
sudo docker compose down
sudo docker compose up -d
```

## Run Tracer
```bash
tracer init --pipeline-name airflow_test --environment sandbox_test --user-operator vincent --pipeline-type rnaseq
```

## How to find the results of a specific pipeline run

1. Unpause the pipeline that you want to run
```bash
sudo docker compose run --rm airflow-cli dags unpause pipeline_rnaseq
```

2. Run that particular pipeline
```bash
sudo docker compose run --rm airflow-cli dags trigger pipeline_rnaseq --run-id=my_custom_run_001
```
3. Get the results of that particular pipeline run
```bash
sudo docker compose run --rm airflow-cli tasks states-for-dag-run pipeline_rnaseq my_custom_run_001
```

## Running Airflow Tasks
```bash
# Reapply changes if you need to trobleshoot
sudo docker compose restart airflow-scheduler
```

```bash
# Unpause the DAG
sudo docker compose run --rm airflow-cli dags unpause pipeline_rnaseq

# Trigger the DAG (will reply false)
sudo docker compose run --rm airflow-cli dags trigger pipeline_rnaseq

# Check the status
sudo docker compose run --rm airflow-cli dags state pipeline_rnaseq

# List runs
sudo docker compose run --rm airflow-cli dags list-runs --dag-id pipeline_rnaseq
```

## Clean up old pipeline runs
```bash
sudo docker compose run --rm airflow-cli airflow dags delete pipeline_rnaseq
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
