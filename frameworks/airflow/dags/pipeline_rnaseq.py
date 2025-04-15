from datetime import timedelta
from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.utils.dates import days_ago

default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=1),
}

with DAG(
    dag_id='pipeline_rnaseq',
    default_args=default_args,
    description='RNA-seq analysis pipeline using external bash scripts',
    schedule_interval=None,
    start_date=days_ago(1),
    catchup=False,
    tags=['rnaseq', 'bioinformatics'],
) as dag:

    script1 = BashOperator(
        task_id='run_script1',
        bash_command='bash "/opt/airflow/scripts/script1.sh"',
    )

    script2 = BashOperator(
        task_id='run_script2',
        bash_command='bash "/opt/airflow/scripts/script2.sh"',
    )

    script3 = BashOperator(
        task_id='run_script3',
        bash_command='bash "/opt/airflow/scripts/script3.sh"',
    )

    script1 >> script2 >> script3
