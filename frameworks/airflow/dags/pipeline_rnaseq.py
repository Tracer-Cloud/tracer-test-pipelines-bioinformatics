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
    dag_id='pipeline_rnaseq',
    default_args=default_args,
    description='One-time RNA-seq analysis pipeline',
    schedule_interval=None,  # Manual trigger only
    start_date=datetime(2025, 1, 1),  # Fixed past date to prevent auto-triggering
    catchup=False,  # Do not backfill missed runs
    tags=['rnaseq', 'manual', 'adhoc'],
) as dag:

    script1 = BashOperator(
        task_id='run_script1',
        bash_command='bash "/opt/airflow/scripts/1.quality_control.sh"',
        do_xcom_push=False,
    )

    script2 = BashOperator(
        task_id='run_script2',
        bash_command='bash "/opt/airflow/scripts/2.alignment.sh"',
        do_xcom_push=False,
    )

    script3 = BashOperator(
        task_id='run_script3',
        bash_command='bash "/opt/airflow/scripts/3.quantification.sh"',
        do_xcom_push=False,
    )

    script1 >> script2 >> script3
