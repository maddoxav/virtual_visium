import os

# not cross-platform!!! windows only right now!!!
# replace with path to your own openslide bin folder
OPENSLIDE_PATH = (
    r"C:\Users\maddoxav\openslide-win64-20171122\openslide-win64-20171122\bin"
)
with os.add_dll_directory(OPENSLIDE_PATH):
    import openslide
    from openslide import open_slide
from pathlib import Path
from PIL import Image
from dataclasses import dataclass
import pandas as pd
import random
import uuid

@dataclass
class TileInfo:
    """Tiling information for a WSI
    Attributes:
        wsi (openslide.OpenSlide): whole slide image
        coords (list[tuple[x, y]]): coordinates for top left pixel of tiles at level 0
        level (int): level in wsi
        size (tuple[width, height]): tile region size
    """

    wsi: openslide.OpenSlide
    coords: list[tuple[int, int]]
    level: int
    size: tuple[int, int]


def wsi_tile_generator(tile_info):
    """Generate tiles from H&E WSI at spatial coordinates
    Args:
        tile_info (TileInfo): information about tiles for a WSI
    Yields:
        tile (PIL.Image): tile image
    """
    for coords in tile_info.coords:
        yield tile_info.wsi.read_region(coords, tile_info.level, tile_info.size).convert("RGB")


def sample_tile_generator(sample_tile_info, shuffle=True, seed=42):
    """Generate tiles from sample of H&E WSIs at spatial coordinates
    Args:
        sample_tile_info (list[TileInfo]): list of tile information for each WSI
        shuffle (bool): shuffle order of tiles
        seed (int): seed for shuffling
    Yields:
        tile (PIL.Image): tile image
    """
    idx_to_wsi = dict()
    sample_tile_info_list = list()
    for i, tile_info in enumerate(sample_tile_info):
        idx_to_wsi[i] = tile_info.wsi
        tile_info_list = [
            (i, coords, tile_info.level, tile_info.size) for coords in tile_info.coords
        ]
        sample_tile_info_list += tile_info_list
    if shuffle:
        random.seed(seed)
        random.shuffle(sample_tile_info_list)
    for (idx, coords, level, size) in sample_tile_info_list:
        yield idx_to_wsi[idx].read_region(coords, level, size).convert("RGB")

def tile_writer(tile_generator, dest_dir):
    """Write tiles from generator to jpg files
    Args:
        tile_generator (generator): generator of tiles
        dest_dir (pathlib.Path): directory where tiles will be saved
    """
    if not dest_dir.exists():
        os.mkdir(dest_dir)

    for tile in tile_generator:
        filename = str(uuid.uuid4().hex) + '.jpg'
        dest_path = dest_dir / filename
        tile.save(dest_path)

def coords_from_file(file_path):
    """Generate list of tile coordinates from file
    Args:
        file_path (pathlib.Path): coordinate file path, xs are 1st col, ys are 2nd col
    Returns:
        coords (list[tuple[x, y]]): coordinates for top left pixel of tiles at level 0
    """
    coord_df = pd.read_csv(file_path, sep=",", header=None)
    return list(coord_df.itertuples(index=False, name=None))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate tiles from wsi file")
    parser.add_argument("--wsi_file", type=str, help="path to wsi file")
    parser.add_argument("--dest_dir", type=str, help="directory where tiles will go")
    parser.add_argument("--coord_file", type=str, help="path to coordinate file")
    parser.add_argument("--level", type=int, help="level of wsi")
    parser.add_argument("--width", type=int, help="tile width")
    parser.add_argument("--height", type=int, help="tile height")
    args = parser.parse_args()

    wsi = open_slide(args.wsi_file)
    coords = coords_from_file(args.coord_file)
    tile_info = TileInfo(wsi, coords, args.level, (args.width, args.height))

    tile_writer(wsi_tile_generator(tile_info), Path(args.dest_dir))
