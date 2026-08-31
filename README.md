
# ggslicer

<!-- badges: start -->

<!-- badges: end -->

This repository is for the ‘ggslicer’ R package.

Scientific visualization is an integral aspect to sharing research data.
In the field of data science, data is commonly visualized with the
‘ggplot2’ package in R. ‘ggslicer’ is a package available for both R and
python that efficiently wrangles imaging data into a data frame to be
utilized by ‘ggplot2’ for image visualization.

Built around the ‘SimpleITK’ image processing framework, ‘ggslicer’ is
able to process 2D and 3D images such as DICOM, Nifti, MINC, and
Analyze.

**How to get started:** There are currently four functions in the
package; two of which may be used to extract a data frame, and two that
fix and preserve MINC image orientation when reading and writing the
image with ‘SimpleITK’. The four functions are of the following:

**slice_axis() :** Read an image file and return a data frame containing
voxel and world coordinates along specified slice.

**slice_intensity() :** Read an image file and optional data frame with
specific coordinates and return a data frame with intensity values
appended to each coordinate.

**ReadImage_fix() :** Read an image file with ‘SimpleITK’.

**WriteImage_fix() :** Write an image file with ‘SimpleITK’.

## Installation
First, ensure you have SimpleITK installed on your system. 

**On Mac:**

CMake is required to compile SimpleITK, which can be done via homebrew:
``` bash
brew install cmake
```
Then, use a remotes based installer for SimpleITK in R with:
``` r
remotes::install_github("SimpleITK/SimpleITKRInstaller")
```
**Note:** More information on how to install SimpleITK in R can be found at SimpleITK/SimpleITKRInstaller (https://github.com/SimpleITK/SimpleITKRInstaller).

Finally, you can install the development version of ‘ggslicer’ in R by simply
running:

``` r
pak::pak("shanafu/ggslicerR")
library(ggslicer)
```

## Example

This is a basic example using the slice_axis() function in ‘ggslicer’:

``` r
library(ggslicer)

file <- "~/path/to/test/file"

coordinates_df <- slice_axis(file, "y", c(7.5, 88))
view(coordinates_df)
```

