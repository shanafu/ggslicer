# ggslicer

This repository is for the `ggslicer` R package.

Scientific visualization is an integral aspect to sharing research data. In the field of data science, data is commonly visualized with the `ggplot2` package in R. `ggslicer` turns slices of medical and scientific images into tidy data frames that `ggplot2` can plot.

Load a DICOM, NIfTI, MINC, or Analyze image, pick a slice (axis-aligned or oblique), and get back a data frame in world coordinates. Segmentations, atlases, contours, grids, and warp fields all land in the same coordinate space, so they overlay directly in one plot.

Built on the `SimpleITK` image processing framework, `ggslicer` handles 2D, 3D, and higher-dimensional (4D/5D) images.

- To find the tutorial, visit: [shanafu/ggslicer-tutorial](https://github.com/shanafu/ggslicer-tutorial)
- To find the equivalent Python package, visit: [shanafu/ggslicer-python](https://github.com/shanafu/ggslicer-python)

<p>
<img width="1344" height="960" alt="Thalamus atlas overlay on a template brain with labeled structures" src="https://github.com/user-attachments/assets/3edfdf7d-1847-4e44-9916-e48b6472276e" />  
<em> An example of a neuroimaging figure created in R entirely using Grammar of Graphics-based code (via <code>ggplot2</code>): an atlas overlay on a template with a specific focus on the thalamus and its labeled structures. </em>
</p>

## Installation

`ggslicer` depends on `SimpleITK`, which is compiled from source. Expect the first build to take 30+ minutes and a few GB of free disk space. You need a C++ toolchain and CMake first.

**On macOS:**

A compiler toolchain and CMake are required to compile SimpleITK, which can be done via the Xcode command line and homebrew, respectively:

```bash
xcode-select --install
brew install cmake
```

**On Linux (Debian/Ubuntu):**

```bash
sudo apt update
sudo apt install build-essential cmake git r-base-dev
```

On Fedora/RHEL, use `sudo dnf install gcc-c++ make cmake git R-devel`.

**On Windows:**

1. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version. It provides the compiler.
2. Install [CMake](https://cmake.org/download/) and choose "Add CMake to the system PATH" in the installer.
3. Install [Git for Windows](https://git-scm.com/download/win).

Restart R after installing so it picks up the PATH changes.

**Then, on any platform**, install `SimpleITK` and `ggslicer` from R:
 
```r
install.packages("remotes")
remotes::install_github("SimpleITK/SimpleITKRInstaller")
remotes::install_github("shanafu/ggslicer")
```

Check that everything loaded:
 
```r
library(SimpleITK)
library(ggslicer)
```

**Note:** More information on how to install `SimpleITK` in R can be found at [SimpleITK/SimpleITKRInstaller](https://github.com/SimpleITK/SimpleITKRInstaller).


## Function reference

Every function below is exported. Use `?function_name` for full documentation.

### Slicing and sampling

| Function | Description |
|---|---|
| `slice_image()` | <ul><li>**The main function.**</li><li> Samples an image (object or file path) along an axis/coordinate or a pre-built geometry. </li><li> Returns a tibble with the image in `value` and any `extra_images` (e.g. a mask) as additional columns. </li><li> Use `interpolator_overrides` to change the default interpolation for any image. </li></ul>|
| `sample_images()` | <ul><li> Samples one or more images onto an existing `SliceGeometry`, `SlicePackage`, or `SlicePackageSet`. </li><li> Supports 4D/5D images through `extra_index`. </li></ul> |
| `build_slice_geometry()` | <ul><li> Builds a `SlicePackageSet` from an image, an out-of-plane `axis` (`1`/`2`/`3`, `"x"`/`"y"`/`"z"`, or `"sagittal"`/`"coronal"`/`"axial"`), and world `coordinate` values. </li></ul> |
| `discrete_data_names()` | <ul><li> Lists names treated as discrete data (`"mask"`, `"label"`, `"atlas"`, `"segmentation"`). </li><li> These are always sampled with nearest-neighbour interpolation. </li></ul> |


### Contours

| Function | Description |
|---|---|
| `slice_contours()` | <ul><li> Intensity contour lines of a continuous image at given `levels`, as vertices ready for `geom_path()`. </li></ul>|
| `slice_label_contours()` | <ul><li> Region outlines of a discrete label image (masks, atlases, segmentations). </li></ul> |
| `suggest_contour_levels()` | <ul><li>  Suggests `levels` by quantile (default) or by density troughs between populations. Not for discrete data. </li></ul> |


### Label annotations

| Function | Description |
|---|---|
| `contour_centroids()` | <ul><li> One label position per connected region (area-weighted centroid) from `slice_label_contours()` output. </li><li> Supports `label_names` and `min_area`. </li></ul> |
| `slice_label_annotations()` | <ul><li> Wrapper that runs `slice_label_contours()` then `contour_centroids()`. </li><li>Use it when you only need the text labels. </li></ul> |
| `slice_label_annotations_layer()` | <ul><li> Returns a ready-made `geom_text()` layer for the output of `slice_label_annotations()`. </li><li>Style with `color`, `size`, `alpha`, `fontface`. </li></ul>|


### Grid overlay

| Function | Description |
|---|---|
| `slice_grid()` | <ul><li> Bounding box and gridlines for a slice, as line vertices. </li><li> Control density with `spacing`, `point_spacing`, and `padding`. </li></ul> |
| `slice_grid_layers()` | <ul><li> Returns a named list of two styled `geom_path()` layers: `box` and `grid`. </li></ul>|


### Transforms and warp field arrows

| Function | Description |
|---|---|
| `read_minc_transform()` | <ul><li>Reads an MNI `.xfm` transform file, including `Grid_Transform` and multi-block files, which `SimpleITK::ReadTransform()` can silently mishandle.</li><li> With `corrected = TRUE` (default) the transform is adjusted to work on points from `ReadImage_fix()`-corrected images. </li><li>Set `corrected = FALSE` for the transform exactly as written in the file. </li></ul>|
| `transform_points()` | <ul><li> Applies a SimpleITK transform to the `x`/`y`/`z` columns of any tidy data frame. </li><li>`invert = TRUE` applies the inverse, which works for affine transforms only.</li></ul>|
| `slice_warp_arrows()` | <ul><li> Displacement arrows from a transform or displacement field on a lattice of points in a slice. </li></ul> |
| `slice_warp_arrows_layer()` | <ul><li>Returns a ready-made `geom_segment()` layer with an arrowhead for the output of `slice_warp_arrows()`. </li><li>Style with `color`, `linewidth`, `alpha`, `arrow_length_inches`, `arrow_type`. </li></ul>|


### I/O

| Function | Description |
|---|---|
| `ReadImage_fix()` | <ul><li> Reads an image with SimpleITK; corrects MINC header orientation so it matches an equivalent NIfTI. </li></ul> |
| `WriteImage_fix()` | <ul><li> Writes an image; undoes the MINC correction when writing back to MINC. </li></ul>|


### Classes (R6)

You usually don't need these unless you want to visualize oblique slices and custom stacks. A bare `SliceGeometry` can be used anywhere a `SlicePackage` is expected.

| Class | Description |
|---|---|
| `SliceGeometry` | <ul><li> A single evenly sampled 2D slice in 3D world space (no intensity data). </li><li> Constructors: `$new()`, `$from_corners()`, `$from_center()`, `$from_normal()`, `$from_bounds()`, `$from_image_axis()`.  </li></ul>|
| `SlicePackage` | <ul><li> A stack of parallel slices along a shared normal. </li><li> Constructors: `$new()`, `$from_extent_k()`, `$from_center_k()`, `$from_slices()`, `$from_image_axis()`. </li></ul>|
| `SlicePackageSet` | <ul><li> A named collection of packages (e.g. several orientations) that produces the final long-format table. </li><li> Constructors: `$new()`, `$from_orthogonal_triplet()`, `$from_slice_packages()`. </li></ul>|
