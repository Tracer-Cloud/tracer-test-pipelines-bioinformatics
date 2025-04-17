# Airflow Testing For Tracer 
- We assume tracer runs in the background
- When testing on EC2 AWS you need to be logged in as an ubuntu user and execute all commands with sudo 

# Instruction with celery we do not need a worker?

## Instructions To Run Airflow With Docker
- You need to use docker
```bash
# ALWAYS DO THIS
 sudo docker-compose down --volumes --remove-orphans
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
# Confirm it is using our own image tracer-bio-airflow:latest
sudo docker ps --format "table {{.Names}}\t{{.Image}}"
```

```bash
# Reapply changes if you need to troubleshoot change the code
sudo docker compose restart
```

```bash
# Check logs of the docker container worker
sudo docker logs -f airflow-airflow-worker-1
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
# getting access to the bash environment inside airflow -- type "exit" to exit
sudo docker exec -u 0 -it airflow-airflow-worker-1 /bin/bash
```
```bash
sudo docker exec -it airflow-airflow-worker-1 python /opt/airflow/dags/pipeline_rnaseq.py
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


## Run Tracer
1. Export aws credentials
```bash
#!/bin/bash
export AWS_ACCESS_KEY_ID=XXXXXXX
export AWS_SECRET_ACCESS_KEY=XXXXXXXX
export AWS_LOG_LEVEL="debug"

echo "AWS environment variables exported."
```
2. Run tracer
```bash
tracer init --pipeline-name airflow_vin --environment sandbox_test --user-operator vincent --pipeline-type rnaseq
```

# Concise Todo List To Get Airflow In GitHub Codespaces

Build Docker Image: docker build -t <local-image> .
Create ECR Repo (if needed): AWS Console/CLI.
Authenticate Docker: aws ecr get-login-password ... | docker login ...
Tag for ECR: docker tag <local-image> <ecr-uri>:<tag>
Push to ECR: docker push <ecr-uri>:<tag>
Create .devcontainer/devcontainer.json: In your GitHub repo.
Set image in devcontainer.json: "image": "<ecr-uri>:<tag>"
Handle ECR Auth in Codespaces (if private):
Public: Skip.
Secrets: Configure AWS_* secrets.
Credential Helper: Advanced setup.
Commit & Push .devcontainer.json.
Create GitHub Codespace: It will pull your ECR image.
Verify environment in Codespace.


# New Approach docker images per pipeline step
### Hisat2
sudo docker pull quay.io/biocontainers/hisat2:2.2.1--h503566f_8

### Star
sudo docker pull quay.io/biocontainers/star:2.7.3a--h5ca1c30_1

### fastqc
sudo docker pull quay.io/biocontainers/fastqc:0.11.7--pl5.22.0_2
