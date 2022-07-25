import feather
from pathlib import Path
import pandas as pd
import numpy as np
from torch.utils.data import Dataset
from torchvision.io import read_image
from PIL import Image
from torchvision.transforms import transforms


def parse_tile_name(tile_name: str):
    """Parse tile name into its components"""
    components = tile_name.split("-")
    sample_id = components[0]
    tissue_id = components[1]
    barcode = components[2] + "." + components[3]
    size = components[4].removesuffix(".jpeg")
    size = size.split("x")
    size = tuple(int(dim) for dim in size)
    return {
        "sample_id": sample_id,
        "tissue_id": tissue_id,
        "barcode": barcode,
        "size": size,
    }


def unparse_tile_name(sample_id: str, tissue_id: str, barcode: str, size: tuple):
    """Create tile name from it's components"""
    barcode = "-".join(barcode.split("."))
    return f"{sample_id}-{tissue_id}-{barcode}-{size[0]}x{size[1]}.jpeg"


def parse_count_mat_fname(count_mat_fname: str):
    """Parse count matrix file name into components"""
    components = count_mat_fname.split("-")
    tissue_id = components[0]
    sample_id = components[1]
    return {"sample_id": sample_id, "tissue_id": tissue_id}


def get_symbol_values(count_mat_dpath: Path, symbol: str, seed: int = 42):
    """
    Get symbol(eg gene) values for all sample+tissue count matrices
    where each count matrix is symbols(first row) x barcodes(colnames[1:])
    """
    count_mat_fpaths = list(count_mat_dpath.glob("*-count_mat.feather"))
    values = pd.Series(dtype="float64")
    for mat_fpath in count_mat_fpaths:
        mat_fname = mat_fpath.name
        name_components = parse_count_mat_fname(mat_fname)
        count_mat = feather.read_dataframe(mat_fpath)
        count_mat.set_index(keys="symbol", inplace=True)
        count_mat = count_mat.loc[symbol, :].T.reset_index()
        count_mat.rename({"index": "barcode", symbol: "value"}, axis=1, inplace=True)
        count_mat["sample_id"] = name_components["sample_id"]
        count_mat["tissue_id"] = name_components["tissue_id"]
        count_mat.set_index(keys=["sample_id", "tissue_id", "barcode"], inplace=True)
        new_values = count_mat["value"]
        values = pd.concat([values, new_values], axis=0)
    return values


class VisiumImageDataset(Dataset):
    def __init__(
        self,
        data_path: Path,
        count_mat_dir: str,
        raw_tile_dir: str,
        tile_size: tuple,
        symbol: str,
        
        # at least ToTensor is needed
        transform=transforms.Compose([transforms.RandomHorizontalFlip(),
                                     transforms.ToTensor(),
                                     transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])]),
        target_transform=None
    ):
        count_mat_dpath = data_path / count_mat_dir
        self.tile_labels = get_symbol_values(count_mat_dpath, symbol)
        self.tile_dpath = data_path / raw_tile_dir
        self.tile_size = tile_size
        self.symbol = symbol
        self.transform = transform
        self.target_transform = target_transform

    def __len__(self):
        return len(self.tile_labels)

    def __getitem__(self, idx):
        sample_id, tissue_id, barcode = self.tile_labels.index[idx]
        tile_name = unparse_tile_name(sample_id, tissue_id, barcode, self.tile_size)
        tile_fpath = self.tile_dpath / tile_name
        image = Image.open(str(tile_fpath))  # RGB, W, H, C
        label = self.tile_labels.iloc[idx].tolist()
        if self.transform:
            image = self.transform(image)
        if self.target_transform:
            label = self.target_transform(label)
        return image, label
