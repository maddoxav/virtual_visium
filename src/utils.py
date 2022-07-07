def scale_coords(coords, lr_size, hr_size):
    """Scale coordinates from low to high resolution image
    Args:
        coords (list[tuple[x, y]]): low res image coordinates
        lr_size (tuple[w, h]): low res image size
        hr_size (tuple[w, h]): high res image size
    Returns:
        (list[tuple[x, y]]): high res image coordinates
    """
    x_scale = hr_size[0] // lr_size[0]
    y_scale = hr_size[1] // lr_size[1]
    return [(x * x_scale, y * y_scale) for (x, y) in coords]

def center_to_corner(coords, patch_size):
    """Change coordinates from spot center to patch corner.
    Args:
        coords (list[tuple[x, y]]): coordinates of Visium spot centers
        size (tuple[w, h]): patch size
    Returns:
        (list[tuple[x, y]]): coordinates of patch corner w/ spot at patch center
    """
    x_shift = size[0] // 2
    y_shift = size[1] // 2
    return [(x - x_shift, y - y_shift) for (x, y) in coords]
