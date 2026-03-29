import numpy as np
import rfcnt

print(f"Using rfcnt from {rfcnt.__file__}")

file = "rfcnt_example.npz"
print(f"Loading data from {file}...")
arr = np.load(file)["arr"]

# input("Press Enter")

cw = np.ptp(arr) / 99
result = rfcnt.rfc(arr, class_width=cw, class_offset=arr.min() - cw/2)
print(result)
