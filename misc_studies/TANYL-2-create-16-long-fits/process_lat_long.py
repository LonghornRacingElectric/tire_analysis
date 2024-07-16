import os
import pandas as pd

def _import_data(path: str) -> list[str, list]:
    local_results = {}
    f = open(path, "r")
    _exclude = ["$", "!"]

    data_entry = False

    for line in f:
        char_0 = line.strip()[0]

        if data_entry and (char_0 not in _exclude):
            line_stripped = line.replace(" ", "")

            if "$" in line_stripped:
                line_stripped = line_stripped[:line_stripped.index("$")]

            line_split = line_stripped.split("=")
            if line_split[1].replace(".", "").replace("-", "").replace("E", "").replace("e", "").replace("+", "").replace("\n", "").isnumeric():
                val = float(line_split[1])
            
            else:
                val = line_split[1]

            local_results[list(local_results.keys())[-1]][line_split[0]] = val

        else:
            if (char_0 in _exclude):
                data_entry = False
                continue
        
            if (char_0 == "["):
                if ("[SHAPE]" in line):
                    continue

                local_results[line.strip()[1:-1]] = {}
                data_entry = True
            
    return local_results

tir_files = os.listdir("./12_psi_tir_files")

all_scalars = []

for file in tir_files:
    current_results = _import_data(f"./12_psi_tir_files/{file}")
    long_coeffs = current_results["LONGITUDINAL_COEFFICIENTS"]
    lat_coeffs = current_results["LATERAL_COEFFICIENTS"]
    
    relevant_keys = ["PC1", "PD1", "PD2", "PD3", "PE1", "PE2", "PE3", "PE4", "PK1", "PK2", "PK3", "PH1", "PH2", "PV1", "PV2", "RB1", "RB2", "RC1", "RE1", "RE2", "RH1"]

    relevant_x_keys = [key[:2] + "X" + key[-1] for key in relevant_keys]
    x_dict = {}

    relevant_y_keys = [key[:2] + "Y" + key[-1] for key in relevant_keys]
    y_dict = {}

    for key in long_coeffs.keys():
        if key in relevant_x_keys:
            x_dict[key] = long_coeffs[key]
    
    for key in lat_coeffs.keys():
        if key in relevant_y_keys:
            y_dict[key] = lat_coeffs[key]
    
    x_vals = list(long_coeffs.values())
    y_vals = list(lat_coeffs.values())

    scalars = []
    
    for i in range(len(relevant_keys)):
        if y_vals[i] != 0:
            scale = x_dict[relevant_keys[i][:2] + "X" + relevant_keys[i][-1]] / y_dict[relevant_keys[i][:2] + "Y" + relevant_keys[i][-1]]
        else:
            scale = "N/A"

        scalars.append(scale)
    
    all_scalars.append(scalars)

output_dict = {}

count = 0
for file in tir_files:
    output_dict[file[:-4] + ".csv"] = all_scalars[count]
    
    count += 1

df = pd.DataFrame(output_dict)

df.to_csv(f"./outputs/scaling_factors.csv")