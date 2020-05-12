# Level - 4 

![version](https://img.shields.io/github/v/release/Geet-George/JOANNE?color=teal&include_prereleases&label=Latest%20JOANNE%20VERSION&style=for-the-badge)


<div style="text-align: justify">

Level-4 is a dataset of circle products.
<!-- 
[File Structure](#file-structure)

[Test Section](#test-section) -->

## File Structure

| **OBJECT**      | **NAME**           | **DESCRIPTION**                                                               | **UNITS**                             | **DIMENSION**    |
| --------------- | ------------------ | ----------------------------------------------------------------------------- | ------------------------------------- | ---------------- |
| **Dimensions**  | `sounding`         | Sonde number                                                                  |                                       | sounding         |
|                 | `circle`           | Circle number                                                                 |                                       | circle           |
|                 | `height`           | Geopotential Height                                                           | m                                     | height           |
| **Coordinates** | `launch_time`      | Time of dropsonde launch                                                      | seconds since 1970-01-01 00:00:00 UTC | circle, sounding |
|                 | `flight_height`    | mean height of the aircraft during the circle                                 | m                                     | circle           |
|                 | `circle_x`         | mean zonal distance from zero longitude of all regressed sondes in circle     | m                                     | circle           |
|                 | `circle_y`         | mean meridional distance from zero latitude of all regressed sondes in circle | m                                     | circle           |
|                 | `circle_radius`    | mean radius of circle                                                         | m                                     | circle           |
|                 | `circle_time`      | mean time of circle                                                           | UTC                                   | circle           |
| **Variables**   | `Platform`         | platform of the flown circle                                                  |                                       | circle           |
|                 | `u`                | intercept value from regressed eastward wind in circle                        | m s-1                                 | circle, height   |
|                 | `dudx`             | zonal gradient of eastward wind                                               | s-1                                   | circle, height   |
|                 | `dudy`             | meridional gradient of eastward wind                                          | s-1                                   | circle, height   |
|                 | `sondes_regressed` | number of sondes regressed                                                    |                                       | circle, height   |
|                 | `v`                | intercept value from regressed northward wind in circle                       | m s-1                                 | circle, height   |
|                 | `dvdx`             | zonal gradient of northward wind                                              | s-1                                   | circle, height   |
|                 | `dvdy`             | meridional gradient of northward wind                                         | s-1                                   | circle, height   |
|                 | `q`                | intercept value from regressed specific humidity in circle                    | kg kg-1                               | circle, height   |
|                 | `dqdx`             | zonal gradient of specific humidity                                           | m-1                                   | circle, height   |
|                 | `dqdy`             | meridional gradient of specific humidity                                      | m-1                                   | circle, height   |
|                 | `T`                | intercept value from regressed temperature in circle                          | degree_Celsius                        | circle, height   |
|                 | `dTdx`             | zonal gradient of temperature                                                 | degree_Celsius m-1                    | circle, height   |
|                 | `dTdy`             | meridional gradient of temperature                                            | degree_Celsius m-1                    | circle, height   |
|                 | `p`                | intercept value from regressed pressure in circle                             | hPa                                   | circle, height   |
|                 | `dpdx`             | zonal gradient of pressure                                                    | hPa m-1                               | circle, height   |
|                 | `dpdy`             | meridional gradient of pressure                                               | hPa m-1                               | circle, height   |
|                 | `D`                | horizontal mass divergence                                                    | s-1                                   | circle, height   |
|                 | `vor`              | horizontal mass vorticity                                                     | s-1                                   | circle, height   |
|                 | `W`                | large-scale atmospheric vertical velocity                                     | m s-1                                 | circle, height   |
|                 | `omega`            | large-scale atmospheric pressure velocity                                     | hPa h-1                               | circle, height   |
|                 | `h_adv_q`          | horizontal advection of specific humidity                                     | kg kg-1 s-1                           | circle, height   |
|                 | `h_adv_T`          | horizontal advection of temperature                                           | degree_Celsius s-1                    | circle, height   |
|                 | `h_adv_p`          | horizontal advection of pressure                                              | hPa s-1                               | circle, height   |
|                 | `h_adv_u`          | horizontal advection of eastward wind                                         | m s-1 s-1                             | circle, height   |
|                 | `h_adv_v`          | horizontal advection of northward wind                                        | m s-1 s-1                             | circle, height   |

<!-- ## Test Section -->
