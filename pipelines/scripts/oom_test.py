import numpy as np
import time

data = []
try:
    while True:
        data.append(np.ones((1024, 1024, 10)))  # ~80MB per loop
        print(f'Memory blocks: {len(data)}')
        time.sleep(0.1)
except MemoryError:
    print('MemoryError caught. OOM simulation successful.') 