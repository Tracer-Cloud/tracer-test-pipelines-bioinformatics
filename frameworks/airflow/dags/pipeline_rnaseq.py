from datetime import datetime, timedelta
from airflow import DAG
from airflow.operators.bash import BashOperator

default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 0,  # Optional: Skip retries for one-time runs
}

with DAG(
    dag_id='pipeline_rnaseq_2',
    default_args=default_args,
    description='One-time RNA-seq analysis pipeline',
    schedule_interval=None,  # Manual trigger only
    start_date=datetime(2025, 1, 1),  # Fixed past date to prevent auto-triggering
    catchup=False,  # Do not backfill missed runs
    tags=['rnaseq', 'manual', 'adhoc'],
) as dag:

    script1 = BashOperator(
        task_id='fastqc_version',
        bash_command='fastqc --version',
        do_xcom_push=False,
    )

    script2 = BashOperator(
        task_id='STAR_version',
        bash_command='STAR --version',
        do_xcom_push=False,
    )

    script3 = BashOperator(
        task_id='hisat2_version',
        bash_command='hisat2 --version',
        do_xcom_push=False,
    )

    script1 >> script2 >> script3
