import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from openslide import OpenSlide
from pathlib import Path
import torch 
import feather
from visium_datasets import parse_tile_name


def get_val_spatial_counts(
    val_data: torch.utils.data.dataset.Subset,
    spatial_counts: pd.DataFrame
):
    """Return spatial counts in validation set"""
    val_tile_labels = val_data.dataset.tile_labels[val_data.indices]
    val_barcodes = [barcode for (_, _, barcode) in val_tile_labels.index]
    wsi_barcodes = list(spatial_counts.index)
    wsi_val_barcodes = list(set(val_barcodes).intersection(wsi_barcodes))
    val_spatial_counts = spatial_counts.loc[wsi_val_barcodes]
    return val_spatial_counts


def plot_wsi_coords(
    wsi: OpenSlide,
    spatial_counts: pd.DataFrame,
    title: str,
    scale_by=40,
):
    """Plot location of coords in wsi"""
    dims = wsi.dimensions
    scaled_dims = tuple(x / scale_by for x in dims)
    img = np.array(wsi.get_thumbnail(scaled_dims))
    x = spatial_counts["imagecol"] / scale_by
    y = spatial_counts["imagerow"] / scale_by
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.scatter(x, y, s=3, c="black")
    plt.title(title)
    return fig, ax


def plot_wsi_values(
    wsi: OpenSlide,
    spatial_counts: pd.DataFrame,
    value_col: str,
    title: str,
    scale_by=40,
):
    """Plot values (eg predicted, ground truth, error) over wsi"""
    dims = wsi.dimensions
    scaled_dims = tuple(x / scale_by for x in dims)
    img = np.array(wsi.get_thumbnail(scaled_dims))
    x = spatial_counts["imagecol"] / scale_by
    y = spatial_counts["imagerow"] / scale_by
    c = spatial_counts[value_col]
    fig, ax = plt.subplots()
    ax.imshow(img)
    sc = ax.scatter(x, y, s=3, c=c)
    plt.title(title)
    fig.colorbar(sc, orientation="vertical")
    return fig, ax


def unparse_tc_name(tissue_id: str, sample_id: str):
    """Tissue coordinate file name from components"""
    return f"{tissue_id}-{sample_id}-tc.feather"


def get_predicted_counts(
    predictions_fpath: Path, data_dpath: Path, sample_id: str, tissue_id: str
):
    """Get predicted counts for a sample's tissue"""
    # tissue coordinates data frame
    tiss_coords_fname = unparse_tc_name(tissue_id, sample_id)
    tiss_coords_fpath = data_dpath / tiss_coords_fname
    tiss_coords = feather.read_dataframe(tiss_coords_fpath)
    tiss_coords["barcode"] = tiss_coords["barcode"].str.replace("-", ".")
    tiss_coords.set_index(keys="barcode", inplace=True)
    # predictions data frame
    predictions = pd.read_csv(predictions_fpath, sep=",")
    comps_list = [parse_tile_name(tile_name) for tile_name in predictions["tile_name"]]
    comps_list = [
        (comp["sample_id"], comp["tissue_id"], comp["barcode"]) for comp in comps_list
    ]
    sample_ids, tissue_ids, barcodes = zip(*comps_list)
    predictions["sample_id"] = sample_ids
    predictions["tissue_id"] = tissue_ids
    predictions["barcode"] = barcodes 
    # merge predictions with tissue coordinates
    predictions.set_index(keys="barcode", inplace=True)
    predictions = predictions[
        predictions["tissue_id"].eq(tissue_id) & predictions["sample_id"].eq(sample_id)
    ]
    predictions = pd.concat(
        [predictions, tiss_coords[["imagerow", "imagecol"]]], axis=1
    )
    return predictions


def plot_train_val_loss(loss_fpath: Path | str, title: str, epochs: int):
    """Plot loss curves for training and validation sets"""
    loss_df = pd.read_csv(loss_fpath, sep="\t")
    train_loss = loss_df["train_loss"]
    valid_loss = loss_df["valid_loss"]
    fig = plt.figure()
    plt.plot(train_loss[:epochs], label="train loss")
    plt.plot(valid_loss[:epochs], label="valid loss")
    plt.legend(loc="upper right")
    plt.title(title)
    plt.show()


def plot_train_val_cor(loss_fpath: Path | str, title: str, epochs: int):
    """Plot correlation for training and validation set"""
    loss_df = pd.read_csv(loss_fpath, sep="\t")
    train_cor = loss_df["train_cor"]
    valid_cor = loss_df["valid_cor"]
    fig = plt.figure()
    plt.plot(train_cor[:epochs], label="train cor")
    plt.plot(valid_cor[:epochs], label="valid cor")
    plt.legend(loc="upper right")
    plt.title(title)
    plt.show()


def get_ground_truth_counts(
    count_mat: pd.DataFrame, tiss_coords: pd.DataFrame, symbol: str
):
    tiss_coords = tiss_coords.copy()
    count_mat = count_mat.copy()
    tiss_coords["barcode"] = tiss_coords["barcode"].str.replace("-", ".")
    tiss_coords.set_index(keys="barcode", inplace=True)
    x = tiss_coords["imagecol"]
    y = tiss_coords["imagerow"]
    count_mat.set_index(keys="symbol", inplace=True)
    c = count_mat.loc[symbol]
    return pd.DataFrame(dict(imagecol=x, imagerow=y, count=c))



