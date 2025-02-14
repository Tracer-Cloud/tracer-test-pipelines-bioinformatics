# Tracer-Nextflow

Nextflow pipelines and infrastructure for Tracer client testing.

## Contents

This repo contains a cloudformation template to set up AWS Batch infrastructure
in the Tracer AWS account as well as submodule links to several pipelines that
can be run either locally via Docker or on AWS Batch.

## Setup

Due to the use of submodules, when cloning this repository include the
`--recurse-submodules` option to ensure that submodules are also initialized and
downloaded.

## CloudFormation

### AWS Batch Environment

The [Batch cloudformation template](./cloudformation/nextflow-batch.yml) sets up
AWS Batch compute environments and job queues for Nextflow pipelines as well as
an S3 Bucket, `s3://tracer-nxf-work`, to serve as a work directory for Nextflow
to store intermediate files. It also creates roles and instance profile for the
Batch jobs to read and write to the Nextflow work bucket and pull from other S3
buckets, such as those containing public test data.

There are two compute environment/queues in this template: one for running
CPU-only jobs, `NextflowCPU`, which uses `c`, `r`, and `m` family intances, and
another for GPU-enabled jobs, `NextflowGPU`, which uses `g4dn` instances.

For more convenient debugging, the cloudformation template also creates an EC2
instance connect endpoint, which can be used to SSH into the Batch compute nodes
without needing an SSH keypair. You can use this endpoint to connect to the
instances via the AWS console or the CLI. Non-admin user will need several IAM
permissions to do so. See this link for more details:
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/connect-using-eice.html.

### Cost Usage Report for Right Size Tests

The [Cost Usage Report cloudformation
template](./cloudformation/cost-usage-report.yml) creates a Cost Usage Report
via AWS Data Export in a new S3 bucket, s3://tracer-cur, which saves the
`tracer-cloud` account AWS resource usage as a parquet table that can be queried
with Athena. This can be used to query the actual cost of using most AWS
resources in high granularity.

It's difficult to get cost estimates for AWS Batch/ECS jobs because they can
share the same instances with various packing behavior. AWS has an option for
Data Export, however, that makes this much easier: [split cost
allocation](https://docs.aws.amazon.com/cur/latest/userguide/split-cost-allocation-data.html),
which assigns the cost to ECS jobs as a fraction of an EC2 instance taking into
account all of the jobs running at the same time with their particular CPU and
memory requirements. The Data Export in the included cloudformation
configuration enables split cost allocation for this purpose.

The template also creates a Glue Crawler for the Cost Usage Report and a named
query that uses the right size test tags to analyze the cost by pipeline and
right size test scenario. The Crawler is configured to run on demand, so after
running the tests it is necessary to manually run it before executing the named
query.

> **N.B.**: Data Exports are refreshed several times a day, but not on demand,
> so it is often necessary to wait up to 12 hours after a run completes before
> the ECS tasks appear in the Data Export bucket.

## Pipelines

- [Sarek](./pipelines/nf-core/sarek/)
  - A pipeline for genomic variant analysis, maintained by nf-core.
  - Includes minimal as well as full-scale test data via the `test` and
    `test_full` profiles, respectively.
  - https://nf-co.re/sarek/
- [Rnaseq](./pipelines/nf-core/rnaseq/)
  - A pipeline for RNA-seq analysis, maintained by nf-core.
  - Also includes minimal as well as full-scale test data via the `test` and
    `test_full` profiles.
  - https://nf-co.re/rnaseq/
- [Proteinfold](./pipelines/nf-core/proteinfold/)
  - A pipeline for protein folding via AlphaFold2 or ESMfold, maintained by
    nf-core.
  - Requires a GPU for running full scale tests. Includes a few basic smoke
    tests for Nextflow functionality via the `test_<tool>` profiles. There are
    also several full scale tests via the `test_full_<tool>` profiles that
    require GPU access and signficant memory.
  - https://nf-co.re/proteinfold/

### Running Pipelines

The pipelines can be invoked using the included [Makefile](Makefile), which has
targets for both local and AWS Batch test execution. The following test targets
are available:

  - `test_sarek`
  - `test_sarek_aws_batch`
  - `test_full_sarek_aws_batch`
  - `test_rnaseq`
  - `test_rnaseq_aws_batch`
  - `test_full_rnaseq_aws_batch`
  - `test_proteinfold`
  - `test_proteinfold_aws_batch`
  - `test_full_proteinfold_aws_batch`

The AWS Batch test targets will launch the jobs via AWS Batch and store the
results locally. Users launching jobs on AWS Batch will require AWS credentials
with access to the Tracer AWS account, as well as IAM permissions for creating
and managing AWS Batch resources.

#### Local Execution

An additional Nextflow configuration,
[nextflow-config/local.config](nextflow-config/local.config), can be used to
limit the number of CPU cores or memory during local execution. It also sets the
platform to `aarch64` when running pipelines with Docker images on an ARM
processor, such as Apple silicon.

The local configuration is included automatically for local test targets in the
Makefile.

## Dependencies

In order to run the pipelines, you must have Java 17+ installed as well as the
`nextflow` executable. To create a reproducible environment for running the
piplines, [Spack](https://spack.io/) is used with a [spack.yaml](./spack.yaml)
that specifies a recent version of Nextflow for testing. Building the
environment automatically pulls in the required Java version and all of its
dependencies.

For running pipelines locally, Docker is also required. On MacOS, [Rancher
Desktop](https://rancherdesktop.io/) is the recommended way to provision the
Docker engine. If using an Arm Mac, you will want to enable x86 emulation via
Rosetta 2 via the Rancher Desktop settings. On a Linux system, Docker should be
installed via the appropriate package manager using the
[instructions](https://docs.docker.com/engine/install/) on the Docker website.

To setup the Spack environment there is a Makefile target, `setup_environment`,
which downloads and initializes Spack in the current directory. It then installs
the required dependencies for running Nextflow. In the Makefile test targets,
the Spack environment is automatically activated before running tests. To
manually activate the Spack environment in your current shell session run:

```sh
$ . spack/share/spack/setup-env.sh
$ spack env activate -d .
```

If the environment is already installed, such as by `make setup_environment`,
this will put `nextflow` on your `PATH`.