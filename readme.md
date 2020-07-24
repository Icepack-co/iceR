
# Installation instructions

This package requires the following packages to already be installed:
* R (>= 3.4.4),
* ggplot2 (>= 3.1.1),
* sf (>= 0.9-0),
* dplyr (>= 0.8.5),
* purrr (>= 0.3.0),
* Rcpp (>= 1.0.3),
* devtools (>= 1.13.6),
* httr (>= 1.4.0),
* RProtoBuf (>= 0.4.12),
* leaflet (>= 2.0.2),
* jsonlite (>= 1.6.1),
* magrittr (>= 1.5),
* scales (>= 1.1.0)

The trickiest package to install here is `sf` because of the dependency on the gdal and other spatial libraries. Check out their github installation instructions [here](https://github.com/r-spatial/sf). We recommend the devtools installation process which requires `devtools` installed in R: `install.packages("devtools")`.

`sf` can then also be installed from github: `install_github("r-spatial/sf")`. The balance of the packages are relatively straight forward to install using `install.packages("packagename")` if the installation doesn't automatically upgrade them for you.

## Windows 

Installation on windows is pretty straight forward if the package dependencies above are met:

```
devtools::install_github("icepack-co/iceR")
```
## Linux

### Prerequisites
You will need `protoc` compiler (version >= 3.6.1) available in your path.

Installing from source is the best way to go about this as package managers typically deliver a version from the 2.x.x series. To compile from source on debian you can run:

* `apt-get install autoconf automake libtool curl make g++ unzip`
* You can download [3.6.1 here](https://github.com/protocolbuffers/protobuf/releases/download/v3.6.1/protobuf-all-3.6.1.tar.gz) or download a newer version if you prefer from the [github releases page](https://github.com/protocolbuffers/protobuf/releases).
* Unzip the contents of protobuf-all-x.x.x
* `./configure`
* `make`
* `make check`
* `sudo make install`
* `sudo ldconfig`


### Package installation

You'll need `devtools` installed in R: `install.packages("devtools")`.

Then you can install the `iceR` package directly from github using: 
```
devtools::install_github("icepack-co/iceR")
```

## Mac

You'll need to install `protoc` (version >= 3.6.1) so that it's available in your path and `devtools` (in R) before installing `iceR`.
After which, you can run the devtools installation:
```
devtools::install_github("icepack-co/iceR")
```

# Using the iceR package

The package assumes there is a `config.json` file which contains the details of the endpoint (typically `https://api.icepack.ai`) and your api-key (available from the [icepack portal](https://portal.icepack.ai)). Once you have a key you can run the examples available in R from the [examples repo](https://github.com/icepack-co/examples).

