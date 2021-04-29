# %%
from joanne.Level_4 import dicts
import joanne

from importlib import reload

reload(dicts)

# %%
dicts.list_of_vars

directory = "/Users/geet/Documents/JOANNE/joanne/Level_4/"

var_name = dicts.list_of_vars

Dimensions = [
    "circle",
]
Coordinates = [
    "launch_time",
    "alt",
    # "flight_height",
    "circle_lat",
    "circle_lon",
    # "circle_radius",
    "segment_id",
    "sounding",
]
Variables = [
    var
    for var in var_name
    if (var not in Dimensions) & (var not in Coordinates) & (var[0:3] != "se_")
]


def string_table_row(var):

    if "long_name" in dicts.nc_attrs[var]:
        desc = dicts.nc_attrs[var]["long_name"]
    else:
        desc = dicts.nc_attrs[var]["standard_name"]
    units = dicts.nc_attrs[var]["units"]
    dims_list = dicts.nc_dims[var]
    dims = ", ".join(dims_list)
    str_trow = f" & {var} & {desc} & {units} & {dims} \\\ \\hline \n"

    return str_trow


def rows_for_objects(Object):

    id_ = 0
    for i in var_name:
        if i in eval(Object):
            if id_ == 0:
                file.write(f"{Object}" + string_table_row(i))
                id_ += 1
            else:
                file.write(string_table_row(i))
                id_ += 1


file = open(f"{directory}latex_table_Level_4_v{joanne.__version__}.txt", "w",)

file.write("\\begin{table}[H]\n")
file.write("\\centering\n")
file.write(
    "\\begin{tabular}{p{0.08\\linewidth} p{0.11\\linewidth} p{0.3\\linewidth} p{0.2\\linewidth} p{0.1\\linewidth}}\n"
)
file.write("\\hline\n")
file.write("OBJECT & NAME & DESCRIPTION & UNITS & DIMENSION \\\ \\hline\n")

for Object in ["Dimensions", "Coordinates", "Variables"]:
    rows_for_objects(Object)

file.write("\\end{tabular}\n")

file.write(
    "\caption{Table shows the structure for the Level-4 product, outlining the dimensions, coordinates, variables and their corresponding descriptions, units and dimensions. The ancillary variables giving standard error of the regression estimates are not shown in the table.}\n"
)
file.write("\label{l4}\n")

file.write("\\end{table}")

file.close()
# %%
