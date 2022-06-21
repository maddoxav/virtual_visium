from pathlib import Path
import os

# not cross-platform!!! windows only right now!!!
# replace with path to your own vips bin folder
# note!!! can not import pyvips and openslide together, some DLL path issue??? 
VIPSHOME = "C:\\Users\\maddoxav\\vips-dev-w64-web-8.12.2\\vips-dev-8.12\\bin"
with os.add_dll_directory(VIPSHOME):
    import pyvips


def tile_pyramidize(tif_path, out_path):
    """tile and pyramidize a tiff image to work with openslide
    Args:
        tif_path (Path): path to tiff to be processed
        out_path (Path): path to new tiff will go
    Returns:
        None
    """
    in_img = pyvips.Image.new_from_file(str(tif_path))
    in_img.tiffsave(
        str(out_path),
        compression="none",
        tile=True,
        tile_width=256,
        tile_height=256,
        pyramid=True,
    )


def tile_pyramidize_dir(tif_dir, dest_dir):
    """tile and pyramidize all tiff files in directory
    Args:
        tif_dir (Path): directory with tiffs to be processed
        dest_dir (Path): destination directory for processed tiffs
    Returns:
        None
    """
    # Create destination directory if it doesn't exist
    if not Path.exists(dest_dir):
        os.mkdir(dest_dir)

    # Process all tiffs in tif_dir
    if not Path.is_dir(tif_dir):
        raise ValueError("Oops! tif_dir must be a directory!")
    for tif_path in tif_dir.glob("*.tif"):
        tif_name = tif_path.name
        dest_path = dest_dir / tif_name
        tile_pyramidize(tif_path, dest_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Tile and pyramidize all tifs in directory"
    )
    parser.add_argument(
        "--tif_dir", type=str, help="directory where untiled/unpyramidized tiffs are"
    )
    parser.add_argument(
        "--dest_dir", type=str, help="directory where tiled/pyramidized tiffs will go"
    )
    args = parser.parse_args()

    tif_dir = Path(args.tif_dir)
    dest_dir = Path(args.dest_dir)
    tile_pyramidize_dir(tif_dir, dest_dir)
