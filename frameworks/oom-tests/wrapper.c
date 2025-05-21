
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    const char *python_bin = "/opt/conda/bin/python";
    const char *script_path = "/root/tracer-test-pipelines-bioinformatics/frameworks/oom-tests/usecase_out_of_memory.py";

    char *const args[] = { "python", (char *)script_path, NULL };
    execv(python_bin, args);

    perror("execv failed");
    return 1;
}
