from datetime import datetime
from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.operators.python import PythonOperator
import logging

# Simple default arguments
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'retries': 1,
}

# Basic logging function
def log_message(message):
    logging.info(f"RNA-seq Pipeline: {message}")

# Create DAG
with DAG(
    dag_id='pipeline_rnaseq',
    default_args=default_args,
    description='RNA-seq analysis tools version check',
    schedule=None,  # Manual trigger only
    start_date=datetime(2025, 1, 1),
    catchup=False,
    tags=['rnaseq'],
) as dag:

    # Start message
    start_log = PythonOperator(
        task_id='start_log',
        python_callable=log_message,
        op_kwargs={'message': 'Pipeline started'},
    )

    # Tool version checks
    fastqc_version = BashOperator(
    task_id='fastqc_version',
    bash_command='export PATH="$PATH:/opt/conda/envs/rnaseq/bin" && fastqc --version',
    )

    star_version = BashOperator(
        task_id='star_version',
        bash_command='export PATH="$PATH:/opt/conda/envs/rnaseq/bin" && STAR --version',
    )

    hisat2_version = BashOperator(
        task_id='hisat2_version',
        bash_command='export PATH="$PATH:/opt/conda/envs/rnaseq/bin" && hisat2 --version',
    )
        
    # End message
    end_log = PythonOperator(
        task_id='end_log',
        python_callable=log_message,
        op_kwargs={'message': 'Pipeline completed'},
    )
    
    # Set task dependencies - run tools in parallel
    start_log >> [fastqc_version, star_version, hisat2_version] >> end_log