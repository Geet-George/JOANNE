# %%

import joanne
import glob
from packaging import version

import xarray as xr

# %%
directory = "/Users/geet/Documents/JOANNE/Data/Level_3/"

lv3_filename = glob.glob(f"{directory}EUREC4A_JOANNE*.nc")

vers = [None] * len(lv3_filename)

for n, i in enumerate(lv3_filename):
    vers[n] = version.parse(i)

lv3 = xr.open_dataset(str(max(vers)))

# Dimensions = []
Coordinates = list(lv3.coords)
Variables = list(lv3.data_vars)

var_name = Coordinates + Variables

# %%


def replace_underscore_for_latex(string):
    underscore = "_"
    latex_underscore = "\_"
    return string.replace(underscore, latex_underscore)


def replace_minus_one_with_superscript(string):
    minus_one = "-1"
    return string.replace(minus_one, "$^{-1}$")


def string_table_row(var, lv3=lv3):

    if "description" in lv3[var].attrs:
        desc = lv3[var].attrs["description"]
    else:
        desc = lv3[var].attrs["long_name"]

    # desc = dicts.nc_attrs[var]["description"]
    if var not in ["launch_time", "interpolated_time", "alt_bnds"]:
        units = lv3[var].attrs["units"]
    elif var == "alt_bnds":
        units = "m"
    else:
        units = "seconds since 2020-01-01"
    dims_list = list(lv3[var].dims)
    dims = ", ".join(dims_list)

    str_trow = f" & {replace_underscore_for_latex(var)} & {replace_underscore_for_latex(desc)} & {replace_minus_one_with_superscript(replace_underscore_for_latex(units))} & {replace_underscore_for_latex(dims)} \\\ \n"

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


# %%

save_directory = "/Users/geet/Documents/JOANNE/joanne/Level_3/"

file = open(f"{save_directory}latex_table_Level_3_v{joanne.__version__}.txt", "w",)

file.write(
    "\\begin{longtable}{p{0.08\\linewidth} p{0.12\\linewidth} p{0.3\\linewidth} p{0.21\\linewidth} p{0.12\\linewidth}}\n"
)
# file.write("\\centering\n")
file.write(
    "\caption{Table shows the structure for the Level-3 product, outlining the coordinates, variables and their corresponding descriptions, units and dimensions.}\n"
)
file.write("\label{tab:l3}\n")
file.write("\\\ \n")
file.write("\\toprule \n")
# file.write(
#     "\\begin{tabular}{p{0.08\\linewidth} p{0.12\\linewidth} p{0.3\\linewidth} p{0.21\\linewidth} p{0.12\\linewidth}}\n"
# )
# file.write("\\hline\n")
file.write("OBJECT & NAME & DESCRIPTION & UNITS & DIMENSION \\\ \n")
file.write("\\midrule \n")
file.write("\\endhead \n")

file.write("\\midrule \n")
file.write("\\multicolumn{5}{r}{{Continued on next page}} \\\ \n")
file.write("\\midrule \n")
file.write("\\endfoot \n")

file.write("\\bottomrule \n")
file.write("\\endlastfoot \n")


for Object in ["Coordinates", "Variables"]:
    rows_for_objects(Object)
    if Object != "Variables":
        file.write("\\midrule \n")

# file.write("\\end{tabular}\n")
# file.write("\\bottomrule \n")

file.write("\\end{longtable}")

file.close()
# %%

