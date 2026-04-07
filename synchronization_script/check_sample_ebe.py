import argparse
from pathlib import Path

import awkward as ak
import numpy as np
import pandas as pd
import uproot


COMPARISON_COLUMNS = [
    "event",
    "run",
    "luminosityBlock",
    "n_jets",
    "n_b_jets",
    "jet_1_pt",
    "jet_1_eta",
    "jet_2_pt",
    "jet_2_eta",
    "jet_2_phi",
    "jet_3_pt",
    "jet_3_eta",
    "jet_3_phi",
    "MET_pt",
    "n_electrons",
    "n_muons",
]

DNA_BRANCH_ALIASES = {
    "event": ["event"],
    "run": ["run"],
    "luminosityBlock": ["luminosityBlock"],
    "n_jets": ["n_jets"],
    "n_b_jets": ["n_b_jets"],
    "jet_1_pt": ["jet_1_pt"],
    "jet_1_eta": ["jet_1_eta"],
    "jet_2_pt": ["jet_2_pt"],
    "jet_2_eta": ["jet_2_eta"],
    "jet_2_phi": ["jet_2_phi"],
    "jet_3_pt": ["jet_3_pt"],
    "jet_3_eta": ["jet_3_eta"],
    "jet_3_phi": ["jet_3_phi"],
    "MET_pt": ["MET_pt"],
    "n_electrons": ["n_electrons"],
    "n_muons": ["n_muons"],
}

N2P_BRANCH_ALIASES = {
    "event": ["event"],
    "run": ["run"],
    "luminosityBlock": ["luminosityBlock", "lumiblock"],
    "n_jets": ["n_jets", "njet"],
    "n_b_jets": ["n_b_jets", "nbdfm"],
    "jet_1_pt": ["jet_1_pt", "j1_pt"],
    "jet_1_eta": ["jet_1_eta", "j1_eta"],
    "jet_2_pt": ["jet_2_pt", "j2_pt"],
    "jet_2_eta": ["jet_2_eta", "j2_eta"],
    "jet_2_phi": ["jet_2_phi", "j2_phi"],
    "jet_3_pt": ["jet_3_pt", "j3_pt"],
    "jet_3_eta": ["jet_3_eta", "j3_eta"],
    "jet_3_phi": ["jet_3_phi", "j3_phi"],
    "MET_pt": ["MET_pt", "met"],
    "n_electrons": ["n_electrons", "nel"],
    "n_muons": ["n_muons", "nmu"],
}

DEFAULT_DNA_DEBUG_BRANCHES = [
    "event",
    "run",
    "luminosityBlock",
    "n_jets",
    "n_b_jets",
    "MET_pt",
    "Jet_pt",
    "Jet_pt_nom",
    "Jet_corrected_pt",
    "Jet_rawFactor",
    "Jet_eta",
    "Jet_phi",
    "Jet_jetId",
    "Jet_jetId_11p9",
    "Jet_neHEF",
    "Jet_neEmEF",
    "Photon_pt",
    "Photon_eta",
    "Photon_phi",
    "Photon_jetIdx",
    "Photon_jetIdx_11p9",
    "Electron_pt",
    "Electron_corrected_pt",
    "Electron_eta",
    "Electron_phi",
    "Muon_pt",
    "Muon_corrected_pt",
    "Muon_eta",
    "Muon_phi",
    "SelectedJet_pt",
    "SelectedJet_eta",
    "SelectedJet_phi",
]

DEFAULT_N2P_DEBUG_BRANCHES = [
    "event",
    "run",
    "lumiblock",
    "met",
    "jet_pt",
    "jet_nanopt",
    "jet_eta",
    "jet_phi",
    "jet_id",
    "jet_ne_emef",
    "jet_isgood",
    "jet_isgood_min",
    "jet_isvetomap",
    "jet_isvetoeta",
    "jet_isphoton",
    "jet_islep",
    "jet_genjet_idx",
    "photon_pt",
    "photon_jet1_dr",
    "photon_jet2_dr",
    "el_pt",
    "el_pt_raw",
    "mu_pt",
    "mu_pt_raw",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare DNA and nano2pico skimmed ntuples event-by-event."
    )
    parser.add_argument("--dna-file", required=True, help="DNA skimmed ROOT file")
    parser.add_argument(
        "--n2p-file",
        help="nano2pico skimmed ROOT file containing both categories",
    )
    parser.add_argument(
        "--n2p-ggf-file",
        help="nano2pico ROOT file used for zero_to_one_jet / ggf category",
    )
    parser.add_argument(
        "--n2p-vbf-file",
        help="nano2pico ROOT file used for two_jet / vbf category",
    )
    parser.add_argument("--output", required=True, help="Output text log path")
    parser.add_argument("--dna-two-jet-tree", default="two_jet")
    parser.add_argument("--dna-zero-one-tree", default="zero_to_one_jet")
    parser.add_argument("--n2p-two-jet-tree", default="two_jet")
    parser.add_argument("--n2p-zero-one-tree", default="zero_to_one_jet")
    parser.add_argument(
        "--mode",
        choices=["both", "n2p_to_dna", "dna_to_n2p"],
        default="both",
        help="Which mismatch direction to dump",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=200,
        help="Maximum number of mismatched events per direction",
    )
    parser.add_argument("--dna-debug-file", help="DNA analysis-level ROOT file")
    parser.add_argument("--n2p-debug-file", help="nano2pico analysis-level ROOT file")
    parser.add_argument(
        "--dna-debug-tree",
        default="test",
        help="DNA analysis tree name used for extra debug dump",
    )
    parser.add_argument(
        "--n2p-debug-tree",
        default="tree",
        help="nano2pico analysis tree name used for extra debug dump",
    )
    parser.add_argument(
        "--dna-debug-event-branch",
        default="event",
        help="Event branch name in DNA analysis tree",
    )
    parser.add_argument(
        "--n2p-debug-event-branch",
        default="event",
        help="Event branch name in nano2pico analysis tree",
    )
    parser.add_argument(
        "--dna-debug-branches",
        default=",".join(DEFAULT_DNA_DEBUG_BRANCHES),
        help="Comma-separated DNA analysis branches to dump for mismatch events",
    )
    parser.add_argument(
        "--n2p-debug-branches",
        default=",".join(DEFAULT_N2P_DEBUG_BRANCHES),
        help="Comma-separated nano2pico analysis branches to dump for mismatch events",
    )
    args = parser.parse_args()

    if not args.n2p_file and not (args.n2p_ggf_file and args.n2p_vbf_file):
        parser.error("Provide either --n2p-file or both --n2p-ggf-file and --n2p-vbf-file")

    return args


def normalize_event_id(value):
    if pd.isna(value):
        return None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        rounded = round(float(value))
        if abs(float(value) - rounded) < 1e-6:
            return int(rounded)
        return float(value)
    try:
        value_float = float(value)
        rounded = round(value_float)
        if abs(value_float - rounded) < 1e-6:
            return int(rounded)
        return value_float
    except (TypeError, ValueError):
        return str(value)


def parse_branch_list(raw_value):
    if not raw_value:
        return []
    return [branch.strip() for branch in raw_value.split(",") if branch.strip()]


def resolve_tree_name(root_file, requested_tree_name):
    if requested_tree_name in root_file:
        return requested_tree_name

    cycle_free_keys = list(root_file.keys(cycle=False))
    if requested_tree_name in cycle_free_keys:
        return requested_tree_name

    if len(cycle_free_keys) == 1:
        return cycle_free_keys[0]

    file_path = getattr(root_file, "file_path", "<unknown file>")
    raise KeyError(
        f"Tree '{requested_tree_name}' not found in {file_path}; "
        f"available trees: {', '.join(cycle_free_keys)}"
    )


def load_tree(file_path, tree_name, branch_aliases):
    with uproot.open(file_path) as root_file:
        resolved_tree_name = resolve_tree_name(root_file, tree_name)
        tree = root_file[resolved_tree_name]
        tree_keys = set(tree.keys())

        selected_branches = {}
        for canonical_name, aliases in branch_aliases.items():
            for alias in aliases:
                if alias in tree_keys:
                    selected_branches[canonical_name] = alias
                    break

        available = [column for column in branch_aliases if column in selected_branches]
        missing = [column for column in branch_aliases if column not in selected_branches]

        df = tree.arrays(list(selected_branches.values()), library="pd").rename(
            columns={source: target for target, source in selected_branches.items()}
        )
    if "event" not in df.columns:
        raise KeyError(f"'event' branch missing in {file_path}:{resolved_tree_name}")
    df["event_key"] = df["event"].map(normalize_event_id)
    df = df.drop_duplicates(subset=["event_key"]).set_index("event_key", drop=False)
    return df, available, missing, resolved_tree_name


def _series_or_default(df, column, default_value):
    if column in df.columns:
        return df[column].fillna(default_value)
    return pd.Series(default_value, index=df.index)


def category_like_mask(df, category_name):
    n_jets = _series_or_default(df, "n_jets", 0)
    n_b_jets = _series_or_default(df, "n_b_jets", 0)
    met_pt = _series_or_default(df, "MET_pt", 0.0)
    n_leptons = _series_or_default(df, "n_electrons", 0) + _series_or_default(df, "n_muons", 0)

    if category_name == "two_jet":
        return (n_jets >= 2) & (n_b_jets == 0) & (n_leptons == 2)
    if category_name == "zero_to_one_jet":
        return (n_jets <= 1) & (met_pt < 90.0) & (n_leptons == 2)
    raise ValueError(f"Unsupported category name: {category_name}")


def fmt_value(value):
    if pd.isna(value):
        return "nan"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def fmt_debug_value(value, max_items=12):
    if isinstance(value, list):
        if len(value) > max_items:
            preview = ", ".join(fmt_debug_value(item, max_items=max_items) for item in value[:max_items])
            return f"[{preview}, ...] (len={len(value)})"
        return "[" + ", ".join(fmt_debug_value(item, max_items=max_items) for item in value) + "]"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6f}"
    if value is None:
        return "None"
    return str(value)


def ak_value_to_python(value):
    python_value = ak.to_list(value)
    if isinstance(python_value, list):
        return python_value
    if isinstance(python_value, (np.integer, int)):
        return int(python_value)
    if isinstance(python_value, (np.floating, float)):
        return float(python_value)
    return python_value


def load_debug_events(file_path, tree_name, event_branch, requested_branches, event_ids):
    if not event_ids:
        return {}, [], requested_branches

    event_keys = {normalize_event_id(event_id) for event_id in event_ids}
    found_rows = {}

    with uproot.open(file_path) as root_file:
        if tree_name not in root_file:
            raise KeyError(f"Tree '{tree_name}' not found in {file_path}")
        tree = root_file[tree_name]
        requested = [event_branch] + [branch for branch in requested_branches if branch != event_branch]
        available = [branch for branch in requested if branch in tree.keys()]
        missing = [branch for branch in requested if branch not in tree.keys()]
        if event_branch not in available:
            raise KeyError(f"'{event_branch}' branch missing in {file_path}:{tree_name}")

        for arrays in tree.iterate(filter_name=available, library="ak", step_size="50 MB"):
            events = ak.to_numpy(arrays[event_branch])
            if len(events) == 0:
                continue
            mask = np.array([normalize_event_id(event_value) in event_keys for event_value in events], dtype=bool)
            if not np.any(mask):
                continue
            selected = arrays[mask]
            selected_events = ak.to_numpy(selected[event_branch])
            for idx, event_value in enumerate(selected_events):
                event_key = normalize_event_id(event_value)
                if event_key in found_rows:
                    continue
                row = {}
                for branch in available:
                    row[branch] = ak_value_to_python(selected[branch][idx])
                found_rows[event_key] = row
                if len(found_rows) == len(event_keys):
                    break
            if len(found_rows) == len(event_keys):
                break

    return found_rows, available, missing


def describe_row(label, row, columns):
    values = [f"{column}={fmt_value(row[column])}" for column in columns if column in row.index]
    return f"{label}:\t" + "\t".join(values)


def describe_diff(source_row, target_row, source_to_target_map):
    diffs = []
    for source_branch, target_branch in source_to_target_map.items():
        if source_branch not in source_row.index or target_branch not in target_row.index:
            continue
        source_value = source_row[source_branch]
        target_value = target_row[target_branch]
        if pd.isna(source_value) and pd.isna(target_value):
            continue
        try:
            delta = target_value - source_value
            if abs(delta) > 1e-6:
                diffs.append(f"{target_branch}-{source_branch}={delta:.6f}")
        except TypeError:
            if target_value != source_value:
                diffs.append(
                    f"{target_branch}-{source_branch}="
                    f"{fmt_value(target_value)} vs {fmt_value(source_value)}"
                )
    if not diffs:
        return "Diff:\tNo significant differences."
    return "Diff:\t" + "\t".join(diffs)


def collect_unmatched_events(source_same_df, target_same_df, limit):
    unmatched = sorted(set(source_same_df.index) - set(target_same_df.index))
    if limit > 0:
        unmatched = unmatched[:limit]
    return unmatched


def describe_debug_row(label, row, columns):
    if row is None:
        return [f"{label}:\tEvent not found in debug tree"]
    lines = [f"{label}:"]
    for column in columns:
        if column in row:
            lines.append(f"  {column}={fmt_debug_value(row[column])}")
    return lines


def resolve_n2p_paths(args):
    shared_path = Path(args.n2p_file) if args.n2p_file else None
    n2p_two_path = Path(args.n2p_vbf_file) if args.n2p_vbf_file else shared_path
    n2p_zero_one_path = Path(args.n2p_ggf_file) if args.n2p_ggf_file else shared_path
    return n2p_two_path, n2p_zero_one_path


def write_direction(
    handle,
    header,
    source_df,
    source_same_df,
    target_same_df,
    target_fallback_df,
    source_label,
    target_label,
    source_same_tree,
    target_same_tree,
    target_fallback_tree,
    source_columns,
    target_columns,
    branch_map,
    unmatched_events,
    source_debug_rows=None,
    target_debug_rows=None,
    source_debug_label=None,
    target_debug_label=None,
    source_debug_columns=None,
    target_debug_columns=None,
):
    handle.write("=" * 100 + "\n")
    handle.write(header + "\n")
    handle.write("=" * 100 + "\n\n")

    for event_id in unmatched_events:
        handle.write(
            f"--- Event {event_id} ({source_label}:{source_same_tree} -> "
            f"{target_label}:{target_same_tree}/{target_fallback_tree}) ---\n"
        )

        source_row = source_df.loc[event_id]
        handle.write(describe_row(f"{source_label} ({source_same_tree})", source_row, source_columns) + "\n")

        if event_id in target_same_df.index:
            target_row = target_same_df.loc[event_id]
            handle.write(describe_row(f"{target_label} ({target_same_tree})", target_row, target_columns) + "\n")
            handle.write(describe_diff(source_row, target_row, branch_map) + "\n\n")
        elif event_id in target_fallback_df.index:
            target_row = target_fallback_df.loc[event_id]
            handle.write(describe_row(f"{target_label} ({target_fallback_tree})", target_row, target_columns) + "\n")
            handle.write(describe_diff(source_row, target_row, branch_map) + "\n\n")
        else:
            handle.write(f"{target_label}:\tEvent not found in {target_same_tree} or {target_fallback_tree}\n\n")

        if source_debug_rows is not None and source_debug_label and source_debug_columns:
            source_debug_row = source_debug_rows.get(event_id)
            handle.write("\n".join(describe_debug_row(source_debug_label, source_debug_row, source_debug_columns)) + "\n")
        if target_debug_rows is not None and target_debug_label and target_debug_columns:
            target_debug_row = target_debug_rows.get(event_id)
            handle.write("\n".join(describe_debug_row(target_debug_label, target_debug_row, target_debug_columns)) + "\n")
        if (
            source_debug_rows is not None and source_debug_label and source_debug_columns
        ) or (
            target_debug_rows is not None and target_debug_label and target_debug_columns
        ):
            handle.write("\n")


def main():
    args = parse_args()

    dna_path = Path(args.dna_file)
    n2p_two_path, n2p_zero_one_path = resolve_n2p_paths(args)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dna_branches = COMPARISON_COLUMNS
    n2p_branches = COMPARISON_COLUMNS

    dna_two_df, dna_two_available, dna_two_missing, dna_two_tree_name = load_tree(
        dna_path, args.dna_two_jet_tree, DNA_BRANCH_ALIASES
    )
    dna_zero_one_df, dna_zero_one_available, dna_zero_one_missing, dna_zero_one_tree_name = load_tree(
        dna_path, args.dna_zero_one_tree, DNA_BRANCH_ALIASES
    )
    n2p_two_df, n2p_two_available, n2p_two_missing, n2p_two_tree_name = load_tree(
        n2p_two_path, args.n2p_two_jet_tree, N2P_BRANCH_ALIASES
    )
    n2p_zero_one_df, n2p_zero_one_available, n2p_zero_one_missing, n2p_zero_one_tree_name = load_tree(
        n2p_zero_one_path, args.n2p_zero_one_tree, N2P_BRANCH_ALIASES
    )

    n2p_two_category_like_df = n2p_two_df[category_like_mask(n2p_two_df, "two_jet")]
    n2p_zero_one_category_like_df = n2p_zero_one_df[
        category_like_mask(n2p_zero_one_df, "zero_to_one_jet")
    ]

    n2p_to_dna_events = []
    dna_to_n2p_events = []
    if args.mode in ("both", "n2p_to_dna"):
        n2p_to_dna_events = collect_unmatched_events(n2p_two_df, dna_two_df, args.limit)
    if args.mode in ("both", "dna_to_n2p"):
        dna_to_n2p_events = collect_unmatched_events(dna_two_df, n2p_two_df, args.limit)

    debug_event_ids = sorted(set(n2p_to_dna_events) | set(dna_to_n2p_events))

    dna_debug_rows = None
    dna_debug_available = []
    dna_debug_missing = []
    if args.dna_debug_file:
        dna_debug_rows, dna_debug_available, dna_debug_missing = load_debug_events(
            file_path=Path(args.dna_debug_file),
            tree_name=args.dna_debug_tree,
            event_branch=args.dna_debug_event_branch,
            requested_branches=parse_branch_list(args.dna_debug_branches),
            event_ids=debug_event_ids,
        )

    n2p_debug_rows = None
    n2p_debug_available = []
    n2p_debug_missing = []
    if args.n2p_debug_file:
        n2p_debug_rows, n2p_debug_available, n2p_debug_missing = load_debug_events(
            file_path=Path(args.n2p_debug_file),
            tree_name=args.n2p_debug_tree,
            event_branch=args.n2p_debug_event_branch,
            requested_branches=parse_branch_list(args.n2p_debug_branches),
            event_ids=debug_event_ids,
        )

    shared_columns = [
        column
        for column in COMPARISON_COLUMNS
        if column in n2p_two_available + n2p_zero_one_available
        and column in dna_two_available + dna_zero_one_available
    ]
    branch_map = {column: column for column in shared_columns}

    with output_path.open("w") as handle:
        handle.write(f"DNA file:\t{dna_path}\n")
        if n2p_two_path == n2p_zero_one_path:
            handle.write(f"n2p file:\t{n2p_two_path}\n")
        else:
            handle.write(f"n2p two_jet file:\t{n2p_two_path}\n")
            handle.write(f"n2p zero_to_one_jet file:\t{n2p_zero_one_path}\n")
        handle.write(f"DNA two_jet entries:\t{len(dna_two_df)}\n")
        handle.write(f"DNA zero_to_one_jet entries:\t{len(dna_zero_one_df)}\n")
        handle.write(f"n2p two_jet entries:\t{len(n2p_two_df)}\n")
        handle.write(f"n2p zero_to_one_jet entries:\t{len(n2p_zero_one_df)}\n")
        handle.write(
            f"n2p two_jet entries comparable to DNA baseline:\t{len(n2p_two_category_like_df)}\n"
        )
        handle.write(
            f"n2p zero_to_one_jet entries comparable to DNA baseline:\t{len(n2p_zero_one_category_like_df)}\n"
        )
        handle.write(
            f"n2p two_jet extra raw-only entries:\t{len(n2p_two_df) - len(n2p_two_category_like_df)}\n"
        )
        handle.write(
            f"n2p zero_to_one_jet extra raw-only entries:\t{len(n2p_zero_one_df) - len(n2p_zero_one_category_like_df)}\n"
        )
        handle.write(
            f"DNA - comparable n2p two_jet:\t{len(dna_two_df) - len(n2p_two_category_like_df)}\n"
        )
        handle.write(
            f"DNA - comparable n2p zero_to_one_jet:\t{len(dna_zero_one_df) - len(n2p_zero_one_category_like_df)}\n"
        )
        handle.write(f"DNA missing branches(two_jet):\t{', '.join(dna_two_missing) or 'none'}\n")
        handle.write(f"DNA missing branches(zero_to_one_jet):\t{', '.join(dna_zero_one_missing) or 'none'}\n")
        handle.write(f"n2p missing branches(two_jet):\t{', '.join(n2p_two_missing) or 'none'}\n")
        handle.write(f"n2p missing branches(zero_to_one_jet):\t{', '.join(n2p_zero_one_missing) or 'none'}\n\n")
        if args.dna_debug_file:
            handle.write(f"DNA debug file:\t{args.dna_debug_file}\n")
            handle.write(f"DNA debug tree:\t{args.dna_debug_tree}\n")
            handle.write(f"DNA debug missing branches:\t{', '.join(dna_debug_missing) or 'none'}\n")
            handle.write(f"DNA debug loaded events:\t{len(dna_debug_rows)}\n")
        if args.n2p_debug_file:
            handle.write(f"n2p debug file:\t{args.n2p_debug_file}\n")
            handle.write(f"n2p debug tree:\t{args.n2p_debug_tree}\n")
            handle.write(f"n2p debug missing branches:\t{', '.join(n2p_debug_missing) or 'none'}\n")
            handle.write(f"n2p debug loaded events:\t{len(n2p_debug_rows)}\n")
        if args.dna_debug_file or args.n2p_debug_file:
            handle.write("\n")

        if args.mode in ("both", "n2p_to_dna"):
            write_direction(
                handle=handle,
                header="Events in n2p two_jet but not in DNA two_jet",
                source_df=n2p_two_df,
                source_same_df=n2p_two_df,
                target_same_df=dna_two_df,
                target_fallback_df=dna_zero_one_df,
                source_label="n2p",
                target_label="DNA",
                source_same_tree=n2p_two_tree_name,
                target_same_tree=dna_two_tree_name,
                target_fallback_tree=dna_zero_one_tree_name,
                source_columns=n2p_branches,
                target_columns=dna_branches,
                branch_map=branch_map,
                unmatched_events=n2p_to_dna_events,
                source_debug_rows=n2p_debug_rows,
                target_debug_rows=dna_debug_rows,
                source_debug_label=f"n2p debug ({args.n2p_debug_tree})" if args.n2p_debug_file else None,
                target_debug_label=f"DNA debug ({args.dna_debug_tree})" if args.dna_debug_file else None,
                source_debug_columns=n2p_debug_available,
                target_debug_columns=dna_debug_available,
            )

        if args.mode in ("both", "dna_to_n2p"):
            write_direction(
                handle=handle,
                header="Events in DNA two_jet but not in n2p two_jet",
                source_df=dna_two_df,
                source_same_df=dna_two_df,
                target_same_df=n2p_two_df,
                target_fallback_df=n2p_zero_one_df,
                source_label="DNA",
                target_label="n2p",
                source_same_tree=dna_two_tree_name,
                target_same_tree=n2p_two_tree_name,
                target_fallback_tree=n2p_zero_one_tree_name,
                source_columns=dna_branches,
                target_columns=n2p_branches,
                branch_map=branch_map,
                unmatched_events=dna_to_n2p_events,
                source_debug_rows=dna_debug_rows,
                target_debug_rows=n2p_debug_rows,
                source_debug_label=f"DNA debug ({args.dna_debug_tree})" if args.dna_debug_file else None,
                target_debug_label=f"n2p debug ({args.n2p_debug_tree})" if args.n2p_debug_file else None,
                source_debug_columns=dna_debug_available,
                target_debug_columns=n2p_debug_available,
            )

    print(f"Wrote mismatch report to {output_path}")


if __name__ == "__main__":
    main()
