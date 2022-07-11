# not cross-platform!!! windows only right now!!!
# replace with path to your own openslide bin folder
import os

OPENSLIDE_PATH = (
    r"C:\Users\maddoxav\openslide-win64-20171122\openslide-win64-20171122\bin"
)
with os.add_dll_directory(OPENSLIDE_PATH):
    from openslide import open_slide
from pathlib import Path
from dataclasses import dataclass
import pandas as pd


def barcode_coords_to_tile_coords(bar_coords, tile_size):
    """Barcode coordinates are visium spot location in low res image from Seurat object.
        Tile coords are top left tile coordinates in high res image with spot at center
        of tile.
    Args:
        bar_coords (pd.DataFrame): barcode, imagerow, imagecol
        tile_size (tuple[w, h]): width and height of tiles
    Returns:
        tile_coords (pd.DataFrame): barcode, imagerow, imagecol
    """
    # shift barcode coord to top left tile coord
    x_shift = tile_size[0] // 2
    y_shift = tile_size[1] // 2
    # modify and return tile coordinates
    tile_coords = bar_coords.copy()
    tile_coords["imagerow"] = bar_coords["imagerow"] - y_shift
    tile_coords["imagecol"] = bar_coords["imagecol"] - x_shift
    return tile_coords


def save_tile(img, dir, sample_id, tissue_id, barcode, tile_size):
    """Save tile as jpeg"""
    tile_size_str = f"{tile_size[0]}x{tile_size[1]}"
    fname = f"{sample_id}-{tissue_id}-{barcode}-{tile_size_str}.jpeg"
    fpath = Path(dir) / fname
    img.save(fpath)


def main(sample_path, wsi_fname, sample_id, tissue_id, tile_size):
    """Tile WSI and save tiles"""
    # read barcode coords df
    barcode_df_fname = f"{tissue_id}-{sample_id}-tc.feather"
    barcode_df_path = Path(sample_path) / "assay_data" / barcode_df_fname
    barcode_df = pd.read_feather(barcode_df_path)
    # read wsi
    wsi_path = Path(sample_path) / "TP_Images" / wsi_fname
    wsi = open_slide(str(wsi_path))
    # tile wsi
    tile_coords = barcode_coords_to_tile_coords(barcode_df, tile_size)
    for _, tc in tile_coords.iterrows():
        row, col = tc['imagerow'], tc['imagecol']
        tile_img = wsi.read_region((col, row), 0, tile_size)
        tile_img = tile_img.convert("RGB")
        tile_dir = Path(sample_path) / "Raw_Tiles"
        save_tile(tile_img, tile_dir, sample_id, tissue_id, tc['barcode'], tile_size)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate tiles from wsi file")
    parser.add_argument("--sample_path", type=str, help="Path to sample directory")
    parser.add_argument("--wsi_fname", type=str, help="WSI file name")
    parser.add_argument("--sample_id", type=str, help="Sample that contains tissue")
    parser.add_argument("--tissue_id", type=str, help="Tissue being tiled")
    parser.add_argument("--tile_width", type=int)
    parser.add_argument("--tile_height", type=int)
    args = parser.parse_args()

    tile_size = (args.tile_width, args.tile_height)

    main(
        args.sample_path,
        args.wsi_fname,
        args.sample_id,
        args.tissue_id,
        tile_size,
    )
