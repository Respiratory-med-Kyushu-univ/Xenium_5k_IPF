import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scFates as scf
import networkx as nx
from scipy.stats import spearmanr
from sklearn.metrics import pairwise_distances

# def hellinger(p, q):
#     p = np.asarray(p, dtype=float)
#     q = np.asarray(q, dtype=float)

#     p = p / p.sum()
#     q = q / q.sum()

#     return np.sqrt(
#         np.sum((np.sqrt(p) - np.sqrt(q))**2)
#     ) / np.sqrt(2)


def held_out_preparation(
    held_out_sample,
    sub_adata,
    score_df,
):
    if isinstance(held_out_sample, str):
        held_out_samples = [held_out_sample]
    else:
        held_out_samples = list(held_out_sample)

    # Exclude held-out samples
    train_mask = ~sub_adata.obs["sample"].isin(held_out_samples)
    train_cells = sub_adata.obs_names[train_mask]

    train_score_df = score_df.loc[train_cells].copy()

    train_adata = anndata.AnnData(
        X=train_score_df.values.copy(),
        obs=sub_adata.obs.loc[train_cells].copy(),
        obsm={
            "spatial": sub_adata[train_cells].obsm["spatial"].copy()
        },
        var=pd.DataFrame(index=train_score_df.columns),
    )

    # Preserve raw neighborhood counts
    train_adata.layers["counts"] = train_adata.X.copy()

    # Preprocessing using training samples only
    sc.pp.filter_cells(train_adata, min_counts=5)
    sc.pp.normalize_total(train_adata, target_sum=1e4)
    sc.pp.log1p(train_adata)

    sc.pp.pca(train_adata, n_comps=5)
    
    # No reference-guided annotation or class exclusion
    return train_adata

    

def trajectory_auto_start(
    traj_adata,
    alveolar_subtypes=(
        "AT1",
        "AT2",
        "Aerocyte",
        "General capillary",
        "Alveolar Pericyte",
    ),
):
    scf.tl.curve(
        traj_adata,
        Nodes=30,
        use_rep="X_pca",
        ndims_rep=2,
        epg_lambda=0.1,
        epg_mu=1,
    )

    cell_xy = traj_adata.obsm["X_pca"][:, :2]
    node_xy = traj_adata.uns["graph"]["F"].T

    dist = np.linalg.norm(
        cell_xy[:, None, :] - node_xy[None, :, :],
        axis=2,
    )

    nearest_node = dist.argmin(axis=1)

    B = traj_adata.uns["graph"]["B"]
    G = nx.from_numpy_array(B)

    tips = traj_adata.uns["graph"]["tips"]

    # Neighborhood proportions using raw counts
    counts = pd.DataFrame(
        traj_adata.layers["counts"],
        index=traj_adata.obs_names,
        columns=traj_adata.var_names,
    )
    
    missing = set(alveolar_subtypes) - set(counts.columns)
    if missing:
        raise ValueError(f"Missing subtypes: {missing}")
    
    alveolar_fraction = (
        counts[list(alveolar_subtypes)].sum(axis=1)
        / counts.sum(axis=1).replace(0, np.nan)
    )

    tip_summary = []

    for tip in tips:
        nodes_near_tip = list(
            nx.single_source_shortest_path_length(
                G,
                int(tip),
                cutoff=3,
            ).keys()
        )

        mask = np.isin(nearest_node, nodes_near_tip)

        local_fraction = alveolar_fraction.iloc[
            np.flatnonzero(mask)
        ]

        tip_summary.append({
            "tip": int(tip),
            "n_cells": int(mask.sum()),
            "alveolar_fraction": local_fraction.mean(),
        })

    tip_summary = pd.DataFrame(tip_summary)

    if tip_summary["alveolar_fraction"].isna().all():
        raise ValueError("Cannot determine root from alveolar fractions.")

    root_tip = int(
        tip_summary.loc[
            tip_summary["alveolar_fraction"].idxmax(),
            "tip",
        ]
    )

    print("Root tip:", root_tip)
    print(tip_summary)

    scf.tl.root(traj_adata, root=root_tip)

    scf.tl.pseudotime(
        traj_adata,
        n_jobs=20,
        n_map=20,
        seed=1,
    )

    t = traj_adata.obs["t"].to_numpy()

    traj_adata.obs["t"] = (
        (t - np.min(t))
        / (np.max(t) - np.min(t))
    )

    return traj_adata

def held_out_projection(n_bins, traj_adata, score_df, held_out_sample, sub_adata, full_cp):
    subtype_columns = score_df.columns.tolist()

    train_df = score_df.loc[traj_adata.obs_names].copy()
    
    train_df["t"] = traj_adata.obs["t"]
    
    ref_prop, bin_centers = make_binned_reference(
        train_df=train_df,
        subtype_columns=subtype_columns,
        continuum_col="t",
        n_bins=n_bins
    )

    if isinstance(held_out_sample, str):
        held_out_samples = [held_out_sample]
    else:
        held_out_samples = list(held_out_sample)

    print("held_out_samples:", held_out_samples)
    
    hos_cells = sub_adata.obs_names[
        sub_adata.obs["sample"].isin(held_out_samples)
    ]
    
    print("Held-out epithelial cells:", len(hos_cells))
    
    hos_cells = hos_cells.intersection(full_cp.index)
    
    print("Held-out continuum cells:", len(hos_cells))
    print("test samples:")
    print(
        sub_adata.obs.loc[hos_cells, "sample"].value_counts()
    )
    
    test_df = score_df.loc[hos_cells].copy()
    print("test_df shape:", test_df.shape)

    
    # Held-out sample(s)の上皮細胞
    # hos_cells = sub_adata.obs_names[sub_adata.obs["sample"].isin(held_out_samples)]
    
    # # Full analysisでcontinuumの対象だったcellだけ
    # hos_cells = hos_cells.intersection(full_cp.index)
    
    # test_df = score_df.loc[hos_cells].copy()
    
    projection_hos = project_to_reference(
        query_df=test_df,
        ref_prop=ref_prop,
        bin_centers=bin_centers,
        subtype_columns=subtype_columns,
        n_nearest_bins=3
    )

    projection_hos["C_full"] = full_cp.loc[projection_hos.index, "C_position"]
    projection_hos["error"] = projection_hos["predicted_continuum_position"] - projection_hos["C_full"]
    projection_hos["absolute_error"] = projection_hos["error"].abs()

    valid = projection_hos[projection_hos["projection_available"]].dropna(
        subset=["predicted_continuum_position","C_full"]
    )
    
    rho, p = spearmanr(valid["C_full"], valid["predicted_continuum_position"])
    
    mae = valid["absolute_error"].mean()
    median_ae = valid["absolute_error"].median()
    median_hd = valid["min_hellinger_distance"].median()

    print("Sample =", held_out_samples)
    print("n =", len(valid))
    print("Spearman rho =", rho)
    print("P =", p)
    print("MAE =", mae)
    print("Median absolute error =", median_ae)
    print("Median min Hellinger =", median_hd)

    return valid, rho, p, mae


# def project_continuum_by_hellinger(
#     ref_df: pd.DataFrame,
#     public_df: pd.DataFrame,
#     subtype_columns: list[str] | None = None,
#     n_nearest_bins: int = 3,
#     eps: float = 1e-8,
# ) -> tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Public dataの上皮細胞ごとのneighborhood compositionを、
#     reference continuum binへHellinger distanceで投影する。

#     Parameters
#     ----------
#     ref_df
#         行がsubtype、列がcontinuum binのDataFrame。
#         各continuum bin列の合計は1を想定。

#     public_df
#         行が上皮細胞、列が周囲subtypeのカウントであるDataFrame。
#         metadata列を含む場合は subtype_columns を指定する。

#     subtype_columns
#         public_df内でsubtypeカウントとして使用する列名。
#         Noneの場合、ref_df.indexとpublic_df.columnsの共通部分を使用する。

#     n_nearest_bins
#         continuum positionの推定に使用する近傍reference bin数。
#         1なら最近傍binのみ、3程度なら距離加重平均。

#     eps
#         距離加重平均でゼロ除算を避ける小さな値。

#     Returns
#     -------
#     projection_df
#         各外部細胞の推定continuum position、最近傍bin、
#         最小Hellinger distanceなど。

#     distance_df
#         各外部細胞と各reference binのHellinger distance。
#     """

#     ref = ref_df.copy()
#     public = public_df.copy()

#     # bin名を数値として扱えるようにする
#     try:
#         ref.columns = ref.columns.astype(float)
#     except (TypeError, ValueError):
#         raise ValueError(
#             "ref_dfの列名はcontinuum binを表す数値に変換可能である必要があります。"
#         )

#     # 使用する共通subtypeを決定
#     if subtype_columns is None:
#         common_subtypes = [
#             subtype
#             for subtype in ref.index
#             if subtype in public.columns
#         ]
#     else:
#         common_subtypes = [
#             subtype
#             for subtype in subtype_columns
#             if subtype in ref.index and subtype in public.columns
#         ]

#     if len(common_subtypes) == 0:
#         raise ValueError("ref_dfとpublic_dfの間に共通subtypeがありません。")

#     print(f"使用する共通subtype数: {len(common_subtypes)}")

#     missing_in_public = [
#         subtype for subtype in ref.index
#         if subtype not in public.columns
#     ]
#     if missing_in_public:
#         print(
#             f"public_dfに存在しないreference subtype数: "
#             f"{len(missing_in_public)}"
#         )

#     # ------------------------------------------------------------
#     # Reference matrix
#     # rows: continuum bins
#     # columns: subtypes
#     # ------------------------------------------------------------
#     ref_matrix = (
#         ref.loc[common_subtypes]
#         .T
#         .sort_index()
#         .astype(float)
#     )

#     # 共通subtypeのみ取り出した後、各binで再度合計1に正規化
#     ref_row_sum = ref_matrix.sum(axis=1)

#     if (ref_row_sum <= 0).any():
#         bad_bins = ref_row_sum.index[ref_row_sum <= 0].tolist()
#         raise ValueError(
#             f"共通subtypeの合計が0となるreference binがあります: {bad_bins}"
#         )

#     ref_prop = ref_matrix.div(ref_row_sum, axis=0)

#     # ------------------------------------------------------------
#     # Public matrix
#     # rows: epithelial cells
#     # columns: subtypes
#     # ------------------------------------------------------------
#     public_counts = public.loc[:, common_subtypes].astype(float)

#     if (public_counts < 0).any().any():
#         raise ValueError("public_dfのsubtype countに負の値があります。")

#     total_neighbor_count = public_counts.sum(axis=1)

#     # 周囲細胞数が0の細胞は投影できない
#     valid_mask = total_neighbor_count > 0

#     public_prop = public_counts.loc[valid_mask].div(
#         total_neighbor_count.loc[valid_mask],
#         axis=0,
#     )

#     # ------------------------------------------------------------
#     # Hellinger distance
#     #
#     # H(p, q) =
#     # 1 / sqrt(2) * ||sqrt(p) - sqrt(q)||_2
#     # ------------------------------------------------------------
#     public_sqrt = np.sqrt(public_prop.to_numpy())
#     ref_sqrt = np.sqrt(ref_prop.to_numpy())

#     distance_array = pairwise_distances(
#         public_sqrt,
#         ref_sqrt,
#         metric="euclidean",
#     ) / np.sqrt(2)

#     distance_df = pd.DataFrame(
#         distance_array,
#         index=public_prop.index,
#         columns=ref_prop.index,
#     )

#     # ------------------------------------------------------------
#     # 最近傍bin
#     # ------------------------------------------------------------
#     nearest_order = np.argsort(distance_array, axis=1)

#     nearest_idx = nearest_order[:, 0]
#     nearest_bin = ref_prop.index.to_numpy()[nearest_idx]
#     min_distance = distance_array[
#         np.arange(distance_array.shape[0]),
#         nearest_idx,
#     ]

#     # 2番目に近いbinとの距離
#     if distance_array.shape[1] >= 2:
#         second_idx = nearest_order[:, 1]
#         second_distance = distance_array[
#             np.arange(distance_array.shape[0]),
#             second_idx,
#         ]
#     else:
#         second_distance = np.full(distance_array.shape[0], np.nan)

#     # ------------------------------------------------------------
#     # 近い複数binによる距離加重平均
#     # ------------------------------------------------------------
#     n_nearest_bins = min(n_nearest_bins, ref_prop.shape[0])

#     selected_idx = nearest_order[:, :n_nearest_bins]

#     selected_distances = np.take_along_axis(
#         distance_array,
#         selected_idx,
#         axis=1,
#     )

#     selected_bins = ref_prop.index.to_numpy()[selected_idx]

#     weights = 1.0 / (selected_distances + eps)
#     weights = weights / weights.sum(axis=1, keepdims=True)

#     predicted_bin = np.sum(
#         weights * selected_bins,
#         axis=1,
#     )

#     # bin 0〜29をcontinuum position 0〜1へ変換
#     min_bin_value = ref_prop.index.min()
#     max_bin_value = ref_prop.index.max()

#     if max_bin_value > min_bin_value:
#         predicted_continuum = (
#             predicted_bin - min_bin_value
#         ) / (
#             max_bin_value - min_bin_value
#         )

#         nearest_continuum = (
#             nearest_bin - min_bin_value
#         ) / (
#             max_bin_value - min_bin_value
#         )
#     else:
#         predicted_continuum = np.zeros_like(predicted_bin)
#         nearest_continuum = np.zeros_like(nearest_bin)

#     # 最小距離と2番目の距離との差
#     # 値が大きいほど最近傍binが比較的明確
#     distance_margin = second_distance - min_distance

#     projection_valid = pd.DataFrame(
#         {
#             "total_neighbor_count": total_neighbor_count.loc[valid_mask],
#             "nearest_bin": nearest_bin,
#             "nearest_continuum_position": nearest_continuum,
#             "predicted_bin": predicted_bin,
#             "predicted_continuum_position": predicted_continuum,
#             "min_hellinger_distance": min_distance,
#             "second_min_hellinger_distance": second_distance,
#             "distance_margin": distance_margin,
#         },
#         index=public_prop.index,
#     )

#     # 周囲細胞数0の細胞も元のindexを維持して返す
#     projection_df = pd.DataFrame(index=public.index)

#     projection_df["total_neighbor_count"] = total_neighbor_count

#     result_columns = [
#         "nearest_bin",
#         "nearest_continuum_position",
#         "predicted_bin",
#         "predicted_continuum_position",
#         "min_hellinger_distance",
#         "second_min_hellinger_distance",
#         "distance_margin",
#     ]

#     projection_df.loc[
#         projection_valid.index,
#         result_columns,
#     ] = projection_valid[result_columns]

#     projection_df["projection_available"] = valid_mask

#     return projection_df, distance_df



def make_binned_reference(
    train_df: pd.DataFrame,
    subtype_columns: list[str],
    continuum_col: str = "continuum_position",
    n_bins: int = 30,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    学習用sampleから、continuumをn_binsに分けたreference compositionを作成する。

    Returns
    -------
    ref_prop
        rows = continuum bins
        columns = subtypes
        各行の合計が1。

    bin_centers
        各binを代表するcontinuum position。
    """

    df = train_df.copy()

    # continuum positionが0〜1の範囲であることを想定
    bin_edges = np.linspace(0, 1, n_bins + 1)

    # 1.0が最後のbinに確実に入るようinclude_lowest=True
    df["_continuum_bin"] = pd.cut(
        df[continuum_col],
        bins=bin_edges,
        labels=False,
        include_lowest=True,
        right=True,
    )

    # 数値誤差などでbinが欠損した細胞を除外
    df = df.dropna(subset=["_continuum_bin"]).copy()
    df["_continuum_bin"] = df["_continuum_bin"].astype(int)

    cell_counts = df[subtype_columns].astype(float)
    cell_total = cell_counts.sum(axis=1)
    
    valid = cell_total > 0
    
    cell_prop = cell_counts.loc[valid].div(
        cell_total.loc[valid],
        axis=0,
    )
    
    cell_prop["_continuum_bin"] = df.loc[
        valid,
        "_continuum_bin",
    ]
    
    ref_prop = (
        cell_prop.groupby("_continuum_bin")[subtype_columns]
        .mean()
        .reindex(range(n_bins))
    )
    
    ref_prop = ref_prop.interpolate(
        method="linear",
        axis=0,
        limit_direction="both",
    )
    
    ref_prop = ref_prop.div(ref_prop.sum(axis=1), axis=0)

    # bin中心
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    ref_prop.index.name = "reference_bin"
    ref_prop.index = np.arange(n_bins)

    return ref_prop, bin_centers


def project_to_reference_ecnp(
    query_df,
    ref_prop,
    bin_centers,
    subtype_columns=None,
    n_nearest_bins=3,
    chunk_size=10000
):
    if subtype_columns is None:
        subtype_columns = ref_prop.columns.tolist()

    subtype_columns = [
        c for c in subtype_columns
        if c in query_df.columns and c in ref_prop.columns
    ]

    query = query_df[subtype_columns].astype(float).copy()
    ref = ref_prop[subtype_columns].astype(float).copy()

    # composition normalization
    query = query.div(query.sum(axis=1), axis=0)
    ref = ref.div(ref.sum(axis=1), axis=0)

    valid = query.notna().all(axis=1)
    query_valid = query.loc[valid]

    Q = np.sqrt(query_valid.values)
    R = np.sqrt(ref.values)

    centers = bin_centers.loc[ref.index].values

    projected = np.full(len(query_valid), np.nan)
    min_distance = np.full(len(query_valid), np.nan)

    for start in range(0, len(query_valid), chunk_size):
        end = min(start + chunk_size, len(query_valid))

        q = Q[start:end]

        # shape:
        # cells × reference bins × cell types
        d = np.sqrt(
            np.sum((q[:, None, :] - R[None, :, :]) ** 2, axis=2) / 2
        )

        min_distance[start:end] = d.min(axis=1)

        # nearest reference bins
        nearest = np.argpartition(
            d,
            kth=n_nearest_bins - 1,
            axis=1
        )[:, :n_nearest_bins]

        nearest_dist = np.take_along_axis(
            d,
            nearest,
            axis=1
        )

        nearest_centers = centers[nearest]

        # inverse-distance weighting
        weights = 1 / (nearest_dist + 1e-8)

        projected[start:end] = (
            np.sum(weights * nearest_centers, axis=1)
            /
            np.sum(weights, axis=1)
        )

    result = pd.DataFrame(
        {
            "projected_position": projected,
            "min_Hellinger": min_distance
        },
        index=query_valid.index
    )

    return result

    
def project_to_reference(
    query_df: pd.DataFrame,
    ref_prop: pd.DataFrame,
    bin_centers: np.ndarray,
    subtype_columns: list[str],
    n_nearest_bins: int = 3,
    eps: float = 1e-8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    query ECNPをHellinger distanceによりreference continuumへ投影する。
    """

    query_counts = query_df[subtype_columns].astype(float).copy()

    if (query_counts < 0).any().any():
        raise ValueError("subtype countに負の値があります。")

    total_neighbor_count = query_counts.sum(axis=1)
    valid_mask = total_neighbor_count > 0

    query_prop = query_counts.loc[valid_mask].div(
        total_neighbor_count.loc[valid_mask],
        axis=0,
    )

    # Hellinger distance
    distance_array = pairwise_distances(
        np.sqrt(query_prop.to_numpy()),
        np.sqrt(ref_prop[subtype_columns].to_numpy()),
        metric="euclidean",
    ) / np.sqrt(2)

    nearest_order = np.argsort(distance_array, axis=1)

    # 最近傍bin
    nearest_idx = nearest_order[:, 0]
    nearest_bin = ref_prop.index.to_numpy()[nearest_idx]

    min_distance = distance_array[
        np.arange(distance_array.shape[0]),
        nearest_idx,
    ]

    # 2番目の距離
    if distance_array.shape[1] >= 2:
        second_idx = nearest_order[:, 1]
        second_distance = distance_array[
            np.arange(distance_array.shape[0]),
            second_idx,
        ]
    else:
        second_distance = np.full(len(nearest_idx), np.nan)

    # 近いn binsを距離加重平均
    n_nearest_bins = min(n_nearest_bins, len(bin_centers))

    selected_idx = nearest_order[:, :n_nearest_bins]

    selected_distances = np.take_along_axis(
        distance_array,
        selected_idx,
        axis=1,
    )

    selected_positions = bin_centers[selected_idx]

    weights = 1 / (selected_distances + eps)
    weights = weights / weights.sum(axis=1, keepdims=True)

    predicted_position = np.sum(
        weights * selected_positions,
        axis=1,
    )

    valid_result = pd.DataFrame(
        {
            "total_neighbor_count": total_neighbor_count.loc[valid_mask],
            "nearest_bin": nearest_bin,
            "nearest_bin_position": bin_centers[nearest_idx],
            "predicted_continuum_position": predicted_position,
            "min_hellinger_distance": min_distance,
            "second_min_hellinger_distance": second_distance,
            "distance_margin": second_distance - min_distance,
        },
        index=query_prop.index,
    )

    # query_dfの全indexを維持
    result = pd.DataFrame(index=query_df.index)

    result["total_neighbor_count"] = total_neighbor_count
    result["projection_available"] = valid_mask

    result_columns = [
        "nearest_bin",
        "nearest_bin_position",
        "predicted_continuum_position",
        "min_hellinger_distance",
        "second_min_hellinger_distance",
        "distance_margin",
    ]

    result.loc[valid_result.index, result_columns] = (
        valid_result[result_columns]
    )

    return result


def assign_ecnp_class(
    query_df: pd.DataFrame,
    class_reference: pd.DataFrame,
    subtype_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    各query ECNPを最も近いECNP classへ割り当てる。
    """

    query_counts = query_df[subtype_columns].astype(float)

    if (query_counts < 0).any().any():
        raise ValueError("subtype countに負の値があります。")

    total_neighbor_count = query_counts.sum(axis=1)
    valid_mask = total_neighbor_count > 0

    query_prop = query_counts.loc[valid_mask].div(
        total_neighbor_count.loc[valid_mask],
        axis=0,
    )

    ref_prop = class_reference.loc[:, subtype_columns].copy()
    ref_prop = ref_prop.div(ref_prop.sum(axis=1), axis=0)

    distance_array = pairwise_distances(
        np.sqrt(query_prop.to_numpy()),
        np.sqrt(ref_prop.to_numpy()),
        metric="euclidean",
    ) / np.sqrt(2)

    # distance_df = pd.DataFrame(
    #     distance_array,
    #     index=query_prop.index,
    #     columns=ref_prop.index,
    # )

    nearest_order = np.argsort(distance_array, axis=1)

    nearest_idx = nearest_order[:, 0]
    second_idx = nearest_order[:, 1]

    predicted_class = ref_prop.index.to_numpy()[nearest_idx]

    min_distance = distance_array[
        np.arange(distance_array.shape[0]),
        nearest_idx,
    ]

    second_distance = distance_array[
        np.arange(distance_array.shape[0]),
        second_idx,
    ]

    second_class = ref_prop.index.to_numpy()[second_idx]

    valid_result = pd.DataFrame(
        {
            "predicted_ecnp_class": predicted_class,
            "class_min_hellinger_distance": min_distance,
            "second_ecnp_class": second_class,
            "class_second_distance": second_distance,
            "class_distance_margin": (
                second_distance - min_distance
            ),
            "total_neighbor_count": (
                total_neighbor_count.loc[valid_mask]
            ),
        },
        index=query_prop.index,
    )

    result = pd.DataFrame(index=query_df.index)

    result["projection_available"] = valid_mask
    result["total_neighbor_count"] = total_neighbor_count

    result.loc[
        valid_result.index,
        [
            "predicted_ecnp_class",
            "class_min_hellinger_distance",
            "second_ecnp_class",
            "class_second_distance",
            "class_distance_margin",
        ],
    ] = valid_result[
        [
            "predicted_ecnp_class",
            "class_min_hellinger_distance",
            "second_ecnp_class",
            "class_second_distance",
            "class_distance_margin",
        ]
    ]

    return result
