# replace with path to your own openslide bin folder
OPENSLIDE_PATH = r'C:\Users\maddoxav\openslide-win64-20171122\openslide-win64-20171122\bin'
import os
with os.add_dll_directory(OPENSLIDE_PATH):
    import openslide 

def tile_generator(wsi, coords, size):
    """ Generate tiles from H&E WSI centered around spatial coordinates
    Args:
        wsi (OpenSlide obj): whole slide image
        coords (list(tuple(int, int))): list of coordinates for tile centers
        size (tuple(int, int)): size of tile in pixels
    Yields:
        tile (np array): tile image
    """
    for (x, y) in coords:
        yield 

