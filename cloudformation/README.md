# How to update the nextflow cloud formation template
```bash
aws cloudformation validate-template --template-body file://nextflow-batch-resources.yml
```


```bash
aws cloudformation deploy \
    --template-file nextflow-batch-resources.yml \
    --stack-name nextflow-batch-resources \
    --capabilities CAPABILITY_NAMED_IAM
```
