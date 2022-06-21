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


@dataclass
class TileInfo:
    """Tiling information for a WSI
    Attributes:
        wsi (openslide.OpenSlide): whole slide image aka WSI
        name (str): name of WSI file
        coords (list[tuple[x, y]]): coordinates for top left pixel of tiles at level 0
        level (int): level in wsi
        size (tuple[width, height]): tile size
    """

    wsi: openslide.OpenSlide
    name: str
    coords: list[tuple[int, int]]
    level: int
    size: tuple[int, int]


@dataclass
class Tile:
    """Tile image and information
    Attributes:
        im (PIL.Image): tile image
        name (str): name of tile's WSI file
        coord (tuple[x,y]): coordinate of tile in wsi
        level (int): level of tile in WSI
        size (tuple[width, height]): tile size
    """

    im: Image
    name: str
    coords: tuple[int, int]
    level: int
    size: tuple[int, int]


def wsi_tile_generator(tile_info):
    """Generate tiles from H&E WSI at spatial coordinates
    Args:
        tile_info (TileInfo): information about tiles for a WSI
    Yields:
        tile (Tile): tile image and info
    """
    for coords in tile_info.coords:
        tile_im = tile_info.wsi.read_region(
            coords, tile_info.level, tile_info.size
        ).convert("RGB")
        tile = Tile(
            im=tile_im,
            name=tile_info.name,
            coords=coords,
            level=tile_info.level,
            size=tile_info.size,
        )
        yield tile


def sample_tile_generator(sample_tile_info, shuffle=True, seed=42):
    """Generate tiles from sample of H&E WSIs at spatial coordinates
    Args:
        sample_tile_info (list[TileInfo]): list of tile information for each WSI
        shuffle (bool): shuffle order of tiles
        seed (int): seed for shuffling
    Yields:
        tile (Tile): tile image and info
    """
    name_to_wsi = dict()
    sample_tile_info_list = list()
    for tile_info in sample_tile_info:
        name_to_wsi[tile_info.name] = tile_info.wsi
        tile_info_list = [
            Tile(None, tile_info.name, coords, tile_info.level, tile_info.size)
            for coords in tile_info.coords
        ]
        sample_tile_info_list += tile_info_list
    if shuffle:
        random.seed(seed)
        random.shuffle(sample_tile_info_list)
    for tile in sample_tile_info_list:
        tile.im = (
            name_to_wsi[tile.name]
            .read_region(tile.coords, tile.level, tile.size)
            .convert("RGB")
        )
        yield tile


def tile_writer(tile_generator, dest_dir):
    """Write tiles from generator to jpg files
    Args:
        tile_generator (generator): generator of tiles
        dest_dir (pathlib.Path): directory where tiles will be saved
    """
    if not dest_dir.exists():
        os.mkdir(dest_dir)

    for tile in tile_generator:
        size_str = f"{tile.size[0]}-{tile.size[1]}"
        coord_str = f"{tile.coords[0]}-{tile.coords[1]}"
        tile_info_str = f"TILE-{coord_str}-{tile.level}-{size_str}"
        filename = f"{tile.name}_{tile_info_str}.jpg"
        dest_path = dest_dir / filename
        tile.im.save(dest_path)


def coords_from_file(file_path):
    """Generate list of tile coordinates from file
    Args:
        file_path (pathlib.Path): coordinate file path, xs are 1st col, ys are 2nd col
    Returns:
        coords (list[tuple[x, y]]): coordinates for top left pixel of tiles at level 0
    """
    coord_df = pd.read_csv(file_path, sep=",", header=None)
    return list(coord_df.itertuples(index=False, name=None))


data_path = r"C:\Users\maddoxav\OneDrive - Michigan Medicine\Documents\Rao_Lab\data\GBM_spatial_transcriptomics\Sample-1512\TP_Images"
tile1 = TileInfo(
    wsi=open_slide(data_path + r"\20H214_Al-Holou_GEX_Slide_Square_A.tif"),
    name="20H214_Al-Holou_GEX_Slide_Square_A",
    coords=[(7700, 7700), (7800, 7800), (7900, 7900)],
    level=0,
    size=(200, 200),
)
tile2 = TileInfo(
    wsi=open_slide(data_path + r"\20H214_Al-Holou_GEX_Slide_Square_C.tif"),
    name="20H214_Al-Holou_GEX_Slide_Square_C",
    coords=[(7700, 7700), (7800, 7800), (7900, 7900)],
    level=0,
    size=(200, 200),
)
sample_tile_info = [tile1, tile2]
tile_writer(
    sample_tile_generator(sample_tile_info),
    Path(
        r"C:\Users\maddoxav\OneDrive - Michigan Medicine\Documents\Rao_Lab\data\GBM_spatial_transcriptomics\tile_images"
    ),
)


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
    tile_info = TileInfo(
        Path(wsi, args.wsi_file).stem, coords, args.level, (args.width, args.height)
    )

    tile_writer(wsi_tile_generator(tile_info), Path(args.dest_dir))
