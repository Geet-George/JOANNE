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

dims = ["time"]
coordinates = ["height", "lat", "lon"]


def string_table_row(var):

    desc = dicts.nc_meta[var]["long_name"]
    units = dicts.nc_meta[var]["units"]
    dims = "time"

    str_trow = f"|{var}|{desc}|{units}|{dims}|\n"

    return str_trow


file.write("|OBJECT|NAME|DESCRIPTION|UNITS|DIMENSION|\n")
file.write(f"|---|---|---|---|---|\n")

id_ = 0
for i in var_name:
    if i in dims:
        if id_ == 0:
            file.write(f"|Dimensions" + string_table_row(i))
            id_ += 1
        else:
            file.write("|" + string_table_row(i))
            id_ += 1

id_ = 0
for i in var_name:
    if i in coordinates:
        if id_ == 0:
            file.write(f"|Coordinates" + string_table_row(i))
            id_ += 1
        else:
            file.write("|" + string_table_row(i))
            id_ += 1

id_ = 0
for i in var_name:
    if (i not in coordinates) & (i not in dims) :
        if id_ == 0:
            file.write(f"|Variables" + string_table_row(i))
            id_ += 1
        else:
            file.write("|" + string_table_row(i))
            id_ += 1

# file.write("|Variable|Long Name|Unit|Dimension|\n")
# file.write("|---|---|---|---|\n")

# for i in range(len(var_name)):
#     file.write(f"|{var_name[i]}|{desc[i]}|{units[i]}|{dims[i]}|\n")

file.close()
# %%
