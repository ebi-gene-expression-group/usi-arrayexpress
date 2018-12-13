
import os

from converter.parsing import parse_sdrf

sdrf_files = ["test_data/" + f for f in os.listdir("test_data") if f.endswith("sdrf.txt")]
for filename in sdrf_files:
    samples, extracts, le, assays, raw_data, processed_data, is_microarray = parse_sdrf(filename)
    for s in samples.values():
        print("Sample: ", s)
    for e in extracts.values():
        print("Extract: ", e)
    for les in le.values():
        print("Lableled Extract: ", les)
    for a in assays.values():
        print("Assay: ", a)
    for f_name, f in raw_data.items():
        print("File: ", f_name, f)
    for m in processed_data.values():
        print("Processed: ", m)
