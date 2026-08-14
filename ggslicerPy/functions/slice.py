import numpy
import itertools
import pandas
from tqdm import tqdm
from .readwrite import ReadImage_fix

def slice_axis(file, slice_orientation, slice_orientation_coords):
    
    image = ReadImage_fix(file)
    file_dims = image.GetSize()

    # Creates all possible voxel combinations
    voxel_indices = list(itertools.product(
        range(file_dims[0]),
        range(file_dims[1]),
        range(file_dims[2])
    ))

    voxel_coords = pandas.DataFrame(voxel_indices, columns=['i', 'j', 'k'])

    # Convert voxel coordinates to world corrdinates using the correct image
    world_coords = []
    for idx, row in tqdm(voxel_coords.iterrows(), total=len(voxel_coords)):
        world_point = image.TransformIndexToPhysicalPoint([int(row['i']), int(row['j']), int(row['k'])])
        world_coords.append(world_point)

    # Convert list to 2D numpy array
    world_coords = numpy.array(world_coords)

    # Add world coordinates as columns
    voxel_coords['x'] = world_coords[:, 0]
    voxel_coords['y'] = world_coords[:, 1]
    voxel_coords['z'] = world_coords[:, 2]

    axis_choice = slice_orientation.lower()

    if axis_choice in ["x", "sagittal", "1"]:
        axis_index = 0
        axis_name = "x"

    elif axis_choice in ["y", "coronal", "2"]:
        axis_index = 1
        axis_name = "y"

    elif axis_choice in ["z", "axial", "3"]:
        axis_index = 2
        axis_name = "z"

    else:
        print("Invalid 'slice_orientation':", slice_orientation, ". Please choose 'x', 'sagittal', '1'; 'y', 'coronal', '2'; or 'z', 'axial', '3'.")
        return None


    filter_col_world = axis_name
    
    available_world_coords = voxel_coords[filter_col_world].unique()

    target_slices = [available_world_coords[numpy/abs(available_world_coords - target_coord).argmin()] for target_coord in slice_orientation_coords]

    target_slices = numpy.unique(target_slices)

    # Filter coordinates based on target slices
    tolerance = 1e-6

    # Create boolean mask
    mask = voxel_coords[filter_col_world].apply(
        lambda coord: numpy.any(numpy.abs(coord - target_slices) < tolerance)
    )

    # Filter and sort
    df_slice_axis = voxel_coords[mask].sort_values([filter_col_world, 'j', 'k'])

    # Reorder columns
    df_slice_axis = df_slice_axis[["x", "y", "z", "i", "j", "k"]].reset_index(drop=True)

    return df_slice_axis


def slice_intensity(file, df_slice_axis=None):
    
    image = ReadImage_fix(file)
    arr = numpy.array(image)

    # Extract intensity values
    intensity_values = numpy.zeros(len(df_slice_axis))

    for row_index in range(len(df_slice_axis)):
        i_coord = df_slice_axis['i'].iloc[row_index]
        j_coord = df_slice_axis['j'].iloc[row_index]
        k_coord = df_slice_axis['k'].iloc[row_index]

        # Check bounds
        if (i_coord < arr.shape[0] and j_coord < arr.shape[1] and k_coord < arr.shape [2] and i_coord >= 0 and j_coord >= 0 and k_coord >= 0):
            intensity_values[row_index] = arr[i_coord, j_coord, k_coord]
        else:
            intensity_values[row_index] = numpy.nan

    df_slice_axis['intensity'] = intensity_values

    # Return properly ordered dataframe
    df = df_slice_axis[["intensity", "x", "y", "z", "i", "j", "k"]]

    return df
