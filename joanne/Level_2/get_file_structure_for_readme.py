# %%
from joanne.Level_2 import dicts
import joanne

from importlib import reload

reload(dicts)

# %%
dicts.list_of_vars

directory = "/Users/geet/Documents/JOANNE/joanne/Level_2/"

file = open(f"{directory}file_structure_v{joanne.__version__}.txt", "w",)

var_name = dicts.list_of_vars

desc = [dicts.nc_meta[var]["long_name"] for var in var_name]
units = [dicts.nc_meta[var]["units"] for var in var_name]
dims = ["time" for _ in var_name]

file.write("|Variable|Long Name|Unit|Dimension|\n")
file.write("|---|---|---|---|\n")

for i in range(len(var_name)):
    file.write(f"|{var_name[i]}|{desc[i]}|{units[i]}|{dims[i]}|\n")

file.close()
# %%
