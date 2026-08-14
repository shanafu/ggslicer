import SimpleITK
from pathlib import Path

def orientation_correction(image):
    flip_filter = SimpleITK.FlipImageFilter()

    # Flip along axes 0 and 1 (first two dimensions)
    flip_axes = [True, True, False] # Flip X and Y, not Z
    flip_filter.SetFlipAxes(flip_axes)

    # Pass image through flip_filter
    flipped_image = flip_filter.Execute(image)

    # Copy metadata
    if image.HasMetaDataKey("OriginalFileType"):
        flipped_image.SetMetaData("OriginalFileType", image.GetMetaData("OriginalFileType"))
    return flipped_image 


def ReadImage_fix(file):
    image = SimpleITK.ReadImage(file)

    # Check image type via extension
    extension = Path(file).suffix.lower()
    is_minc = extension in [".mnc", ".minc"]

    # Flip image orientation if image is MINC
    if is_minc:
        image.SetMetaData("OriginalFileType", "MINC")

        # Store original MINC direction matrix
        orig_direction = image.GetDirection()
        image.SetMetaData("OriginalDirection", ",".join(map(str, orig_direction)))

        # Apply orientation correction for proper display
        corrected_image = orientation_correction(image)

        # Copy metadata to corrected image
        corrected_image.SetMetaData("OriginalFileType", "MINC")
        corrected_image.SetMetaData("OriginalDirection", ",".join(map(str, orig_direction)))

        image = corrected_image

        return image

    # If image type is all else, do not flip
    else:
        image.SetMetaData("OriginalFileType", "Other")

        return image


def WriteImage_fix(image, output_file):
    # Check metadata
    original_file_type_meta = image.HasMetaDataKey("OriginalFileType")
    was_original_minc = False

    if original_file_type_meta:
        original_type = image.GetMetaData("OriginalFileType")
        was_original_minc = original_type == "MINC"

    output_extension = Path(output_file).suffix.lower()
    is_output_minc = output_extension in [".mnc", ".minc"]
    is_output_nifti = output_extension in [".nii"]

    image_to_write = image

    if was_original_minc and is_output_minc:
        # MINC -> MINC: Undo the reading correction to restore proper MINC orientation
        image_to_write = orientation_correction(image)

    elif was_original_minc and is_output_nifti:
        # MINC -> non-MINC: Undo the reading correction to restore proper MINC orientation
        image_to_write = SimpleITK.Image(image) # Copy the image
        ras_direction = [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        image_to_write.SetDirection(ras_direction)

    SimpleITK.WriteImage(image_to_write, output_file)

