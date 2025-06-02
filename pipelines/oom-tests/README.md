# Out-of-Memory (OOM) Test Utility

This folder contains a simple utility to trigger an out-of-memory (OOM) condition, useful for validating Tracer's ability to detect and correlate OOM-killed processes.

## Contents

* `usecase_out_of_memory.py`: A Python script that continuously allocates memory in 100MB blocks until it is terminated by the kernel.
* `wrapper.c`: A minimal C program that wraps the Python script. This helps ensure that Tracer's target matching can properly capture the event as the child process.
* `Makefile` target: Simplifies building and running the test.

## How It Works

The Python script gradually consumes memory in a loop:

```python
import time
blocks = []
try:
    while True:
        blocks.append(bytearray(100 * 1024 * 1024))
        print(f"Allocated {len(blocks) * 100} MB")
        time.sleep(0.1)
except MemoryError:
    print("Caught MemoryError")
    time.sleep(10)
```

The C wrapper simply launches this script using `execv`, passing control fully to Python:

```c
const char *python_bin = "/opt/conda/bin/python";
const char *script_path = "/root/tracer-test-pipelines-bioinformatics/integrations/oom-tests/usecase_out_of_memory.py";
char *const args[] = { "python", (char *)script_path, NULL };
execv(python_bin, args);
```

## Usage

From the `integrations/oom-tests/`:

```bash
make test_out_of_memory_test
```

This will:

1. Compile `wrapper.c` to produce the binary `oom_example_c`
2. Run the binary, which executes the Python script and begins memory allocation

## Why the Wrapper?

Using a C wrapper allows Tracer to correctly identify and match the triggering process, especially in cases where Tracer is configured to filter by executable name or command hierarchy.

## Notes

* Ensure `/opt/conda/bin/python` exists and has the necessary Python environment. Adjust the wrapper path if needed.
* It's recommended to run this on an isolated system or container to avoid destabilizing your host.

---

This utility is designed to intentionally cause system pressure. Monitor system behavior and use with caution.
