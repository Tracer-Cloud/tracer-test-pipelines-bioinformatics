import time

# List to hold large memory blocks
blocks = []

try:
    while True:
        # Allocate 100MB at a time
        blocks.append(bytearray(100 * 1024 * 1024))
        print(f"Allocated {len(blocks) * 100} MB")
        time.sleep(0.1)
except MemoryError:
    print("Caught MemoryError")
    time.sleep(10)  # Keep process alive a bit for observation
