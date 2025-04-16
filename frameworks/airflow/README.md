# Airflow Testing For Tracer 
- We assume tracer runs in the background
- When testing on EC2 AWS you need to be logged in as an ubuntu user and execute all commands with sudo 


## Instructions To Run Airflow With Docker
- You need to use docker

```bash
# The airflow init is only the first time
sudo docker-compose up airflow-init 
```

```bash
# It is very important to run all the different airflow services that we need to run Airflow
sudo docker compose up -d 
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
# Reapply changes if you need to troubleshoot change the code
sudo docker compose restart
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

## Troubleshooting in docker
### Troubleshooting the script inside the environment
```bash
sudo docker exec -it airflow-airflow-worker-1 python /opt/airflow/dags/pipeline_rnaseq.py
```
- sudo docker compose exec -it airflow-worker bash





## Clean up remove logs

```bash
rm -rf logs/* 
```

## Webserver is exposed at:
http://127.0.0.1:8080/login/

Password and username:
- Username: airflow
- Password: airflow


## Run Tracer
```bash
tracer init --pipeline-name airflow_test --environment sandbox_test --user-operator vincent --pipeline-type rnaseq
```

