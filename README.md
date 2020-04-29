# JOANNE - The EUREC<sup>4</sup>A Dropsonde Dataset

JOANNE is the dataset for all dropsondes launched as part of the EUREC<sup>4</sup>A-ATOMIC campaign held in January-February, 2020.

The full form of the acronym JOANNE stands for **J**oint dropsonde-**O**bservations of the **A**tmosphere in tropical **N**orth atla**N**tic large-scale **E**nvironments

---

The different products of JOANNE are listed in the table below. Click on a product to go to its directory, which contains sample data files, along with more information about the product.

| Level                           | Description                                                                                    |
| ------------------------------- | ---------------------------------------------------------------------------------------------- |
| [Level - 0](joanne/Level_0/)           | The raw files of all dropsonde launches from EUREC<sup>4</sup>A                                |
| [Level - 1](joanne/Level_1/)           | Files generated from the ASPEN-processing of all the raw files in Level-0                      |
| [Level - 2](joanne/Level_2/)           | Sounding files that passed additional QC tests, and exclude all soundings with no usable data  |
| [Level - 3A](joanne/Level_3/Level_3A/) | All sounding files in Level-2 gridded on a uniform, vertical grid, with some derived variables |
| [Level - 3B](joanne/Level_3/Level_3B/) | Circle products from all circles flown during EUREC<sup>4</sup>A                               |

---
N.B. : Currently, this git repository is not equipped to serve as a lone-standing package. Think of it as scripts storage for the eventual package. Your feedback (raising issues, requesting features, etc.) will help me speed up this process. :)
