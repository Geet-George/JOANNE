import os
import joanne
import glob
from pylab import size
import joanne.Level_2.QC as qc


data_directory = "/Users/geet/Documents/JOANNE/Data/"
save_directory = "/Users/geet/Documents/JOANNE/Data/QC/"
qc_dir = os.path.join(data_directory, "QC")

jo_version = joanne.__version__

halo_status_file = sorted(glob.glob(f"{qc_dir}*HALO*{jo_version}.nc"))
p3_status_file = sorted(glob.glob(f"{qc_dir}*P3*{jo_version}.nc"))

if (len(halo_status_file) >= 1) and (len(p3_status_file) >= 1):
    print("QC Status files of the current version found for both P3 and HALO")
else:
    print(
        "QC Status files of the current version not found for both P3 and HALO. Therefore, running QC..."
    )

    qc.run_qc(
        data_directory=data_directory, save_directory=save_directory,
    )


l2_files = sorted(glob.glob(f"{data_directory}Level_2/*{jo_version}.nc"))

no_of_l2_files = size(l2_files)

if no_of_l2_files > 1000:
    print(
        f"{no_of_l2_files} Level-2 files found with the current JOANNE version. Therefore not running Level-2 script again."
    )
else:
    print("Starting Level-2")
    import joanne.Level_2.Level_2

    print("Level-2 finished")


l3_files = sorted(glob.glob(f"{data_directory}Level_3/*{jo_version}.nc"))

if size(l3_files) > 0:
    print(
        "Level-3 file found with the current JOANNE version. Therefore not running Level-3 script again."
    )
else:
    print("Starting Level-3")
    import joanne.Level_3.Level_3

    print("Level-3 finished")


l4_files = sorted(glob.glob(f"{data_directory}Level_4/*{jo_version}.nc"))

if size(l4_files) > 0:
    print(
        "Level-4 file found with the current JOANNE version. Therefore not running Level-4 script again."
    )
else:
    print("Starting Level-4")
    import joanne.Level_4.Level_4

    print("Level-4 finished")
