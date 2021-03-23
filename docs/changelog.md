# Changelog for JOANNE

## **VERSIONS WITH v0.10**

### version-in-process

- Scripts added to create $\LaTeX$ tables out of the xarray data structure for Levels 2, 3 and 4 *[[commit](https://github.com/Geet-George/JOANNE/commit/d9c960114f67a38047b753ca7d6506a02f7d75d0)]*
- Variables of standard errors added as ancillary variables to Level-4 *[[commit](https://github.com/Geet-George/JOANNE/commit/f7be804f929fe92def8573f75343ddcad9613736)]*
- Functions now included that will add standard error estimates of the regression variables to Level-4 dataset *[[commit](https://github.com/Geet-George/JOANNE/commit/6637e3fd5839f06900fde6ecddfd00a0f64d7022)]*

### v0.10.1

- $\omega$ no longer provided as variable in Level-4. It is redundant keeping in mind $W$ is available, and density can be estimated using other variables available in JOANNE. The relevant code pieces have also been removed from the scripts *[[commit](https://github.com/Geet-George/JOANNE/commit/0694b6136db896b541022126508a0ce0705d1c53)]*
- Temporary fix for the mismatch between FPS and JOANNE files. This is done by selecting only sondes available in Level-3 even if there are more in FPS. Once `sondes.yaml` and consequently all [FPS YAML files](https://github.com/eurec4a/flight-phase-separation/tree/master/flight_phase_files) are in sync between JOANNE and FPS repo, this fix will be removed and the function `get_circles` from the `ready_ds_for_regression` module will go back to as it is *[[commit](https://github.com/Geet-George/JOANNE/commit/ebaaa818d61c7138e1087d5c7c0847093e1dab9b)]*

### v0.10.0

- `alt` type changes from `int64` to `short`, to make the dataset compatible with OpenDAP functionality *[[commit](https://github.com/Geet-George/JOANNE/commit/cf3d1921764cc41194c4d1c71a2b4ad66ecbaa0e)]*
- The `alt` dimension now has coordinates of altitude values, as opposed to  indices along the altitude dimension, which was due to a bug in the previous major version. The dimension now also has the attribute `bounds` to link it to its cell boundaries' variable `alt_bnds` *[[commit](https://github.com/Geet-George/JOANNE/commit/f517fd49c2d7dc903e1a07a9c9ae14035756123a)]*
- `sonde_id` (with unique IDs) is now the dimension of trajectories for Level-3, as opposed to the `sounding` dimension till now, which were indices along the dimension *[[commit](https://github.com/Geet-George/JOANNE/commit/c304e4039963791db3b3b5e80ef21697fb0ad7a0)]*
- sondes.yaml is now included with Level-2, and is used as an input file for creating the [FPS files](https://github.com/eurec4a/flight-phase-separation). The script to create this YAML file from the status file is also now included in the repo. *[[commit](https://github.com/Geet-George/JOANNE/commit/0f028f22931aef99c65f82d1579dcc3a8ea67cfb)]* 

---
## **VERSIONS WITH v0.9**

### v0.9.4

- `AVAPS Software Notes` and `AVAPS Format Notes` are no longer available as global attributes in Level-2 files *[[commit](https://github.com/Geet-George/JOANNE/commit/5a57201fe90cefe22ce65a235bfffe0f5f291ac9)]*
- `sonde_id` is now added to the status file from Level-2. The sequencing of `sonde_id` is also now chronologically based on `launch_time` rather than the time retrieved from filenames of Level-1 files. This changes the sequence of and around sondes that detected a late launch or failed to detect one at all. *[[commit](https://github.com/Geet-George/JOANNE/commit/98b8aa1ff9ec71edaeb21a0933287f542638e786)]*
- The status file from Level-2 now includes `launch_time` as a variable, which also replaces `time` as dimension along which sonde flags are provided. The variable `launch_time` is also now taken from the attribute of Level-1 files, rather than minimum of `time` values from the Level-1 files *[link to commits [here](https://github.com/Geet-George/JOANNE/commit/f02cdb9de50c2149edef9cce9afb0c51b3c5a429) and [here](https://github.com/Geet-George/JOANNE/commit/e9d147f9cf1961887f20ac3810cf9dbae331d7a1)]*

### v0.9.3

- Sondes for circles selected based on `sonde_id`, and not `launch_time` of segments from FPS files *[[commit](https://github.com/Geet-George/JOANNE/commit/b15599feec221e59416eb0ec741ba702850639f8)]*
- The `fit2d` function now also gives as output the individual sounding profiles of a circle. These will be, in later versions, helpful for estimating uncertainties in Level-4, but the individual soundings will not feature in the Level-4 dataset *[link to commits [here](https://github.com/Geet-George/JOANNE/commit/7d969d28ab6bd22efd8e4e02704f9772c5b64f3b) and [here](https://github.com/Geet-George/JOANNE/commit/e75a37c421c13108a64f6483d9afa467e4ceb306)]*

### v0.9.2

For changes in and before this version, refer to comments on tags in the [JOANNE GitHub Repo](https://github.com/Geet-George/JOANNE)

## Abbreviations

|Abbreviation|Full form / Description|
|:---:|:---:|
|FPS|[Flight phase segmentation files](https://github.com/eurec4a/flight-phase-separation)|
|$W$| meso-scale vertical velocity |
|$\omega$| meso-scale pressure velocity |

