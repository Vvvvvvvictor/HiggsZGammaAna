import argparse
from pathlib import Path

import pandas as pd
import uproot


DEFAULT_BRANCH_MAP = {
    "event": "event",
    "run": "run",
    "luminosityBlock": "luminosityBlock",
    "n_jets": "n_jets",
    "nbdfm": "n_b_jets",
    "jet_1_pt": "jet_1_pt",
    "jet_1_eta": "jet_1_eta",
    "jet_2_pt": "jet_2_pt",
    "jet_2_eta": "jet_2_eta",
    "j2_phi": "jet_2_phi",
    "j3_pt": "jet_3_pt",
    "j3_eta": "jet_3_eta",
    "j3_phi": "jet_3_phi",
    "met": "MET_pt",
    "nel": "n_electrons",
    "nmu": "n_muons",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare DNA and nano2pico skimmed ntuples event-by-event."
    )
    parser.add_argument("--dna-file", required=True, help="DNA skimmed ROOT file")
    parser.add_argument("--n2p-file", required=True, help="nano2pico skimmed ROOT file")
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
    return parser.parse_args()


def load_tree(file_path, tree_name, requested_branches):
    with uproot.open(file_path) as root_file:
        if tree_name not in root_file:
            raise KeyError(f"Tree '{tree_name}' not found in {file_path}")
        tree = root_file[tree_name]
        available = [branch for branch in requested_branches if branch in tree.keys()]
        missing = [branch for branch in requested_branches if branch not in tree.keys()]
        df = tree.arrays(available, library="pd")
    if "event" not in df.columns:
        raise KeyError(f"'event' branch missing in {file_path}:{tree_name}")
    df = df.drop_duplicates(subset=["event"]).set_index("event", drop=False)
    return df, available, missing


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
    limit,
):
    handle.write("=" * 100 + "\n")
    handle.write(header + "\n")
    handle.write("=" * 100 + "\n\n")

    unmatched = sorted(set(source_same_df["event"]) - set(target_same_df["event"]))
    if limit > 0:
        unmatched = unmatched[:limit]

    for event_id in unmatched:
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
            continue

        if event_id in target_fallback_df.index:
            target_row = target_fallback_df.loc[event_id]
            handle.write(describe_row(f"{target_label} ({target_fallback_tree})", target_row, target_columns) + "\n")
            handle.write(describe_diff(source_row, target_row, branch_map) + "\n\n")
            continue

        handle.write(f"{target_label}:\tEvent not found in {target_same_tree} or {target_fallback_tree}\n\n")


def main():
    args = parse_args()

    dna_path = Path(args.dna_file)
    n2p_path = Path(args.n2p_file)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dna_branches = list(DEFAULT_BRANCH_MAP.values())
    n2p_branches = list(DEFAULT_BRANCH_MAP.keys())

    dna_two_df, dna_two_available, dna_two_missing = load_tree(
        dna_path, args.dna_two_jet_tree, dna_branches
    )
    dna_zero_one_df, dna_zero_one_available, dna_zero_one_missing = load_tree(
        dna_path, args.dna_zero_one_tree, dna_branches
    )
    n2p_two_df, n2p_two_available, n2p_two_missing = load_tree(
        n2p_path, args.n2p_two_jet_tree, n2p_branches
    )
    n2p_zero_one_df, n2p_zero_one_available, n2p_zero_one_missing = load_tree(
        n2p_path, args.n2p_zero_one_tree, n2p_branches
    )

    branch_map = {
        n2p_branch: dna_branch
        for n2p_branch, dna_branch in DEFAULT_BRANCH_MAP.items()
        if n2p_branch in n2p_two_available + n2p_zero_one_available
        and dna_branch in dna_two_available + dna_zero_one_available
    }

    with output_path.open("w") as handle:
        handle.write(f"DNA file:\t{dna_path}\n")
        handle.write(f"n2p file:\t{n2p_path}\n")
        handle.write(f"DNA two_jet entries:\t{len(dna_two_df)}\n")
        handle.write(f"DNA zero_to_one_jet entries:\t{len(dna_zero_one_df)}\n")
        handle.write(f"n2p two_jet entries:\t{len(n2p_two_df)}\n")
        handle.write(f"n2p zero_to_one_jet entries:\t{len(n2p_zero_one_df)}\n")
        handle.write(f"DNA missing branches(two_jet):\t{', '.join(dna_two_missing) or 'none'}\n")
        handle.write(f"DNA missing branches(zero_to_one_jet):\t{', '.join(dna_zero_one_missing) or 'none'}\n")
        handle.write(f"n2p missing branches(two_jet):\t{', '.join(n2p_two_missing) or 'none'}\n")
        handle.write(f"n2p missing branches(zero_to_one_jet):\t{', '.join(n2p_zero_one_missing) or 'none'}\n\n")

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
                source_same_tree=args.n2p_two_jet_tree,
                target_same_tree=args.dna_two_jet_tree,
                target_fallback_tree=args.dna_zero_one_tree,
                source_columns=n2p_branches,
                target_columns=dna_branches,
                branch_map=branch_map,
                limit=args.limit,
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
                source_same_tree=args.dna_two_jet_tree,
                target_same_tree=args.n2p_two_jet_tree,
                target_fallback_tree=args.n2p_zero_one_tree,
                source_columns=dna_branches,
                target_columns=n2p_branches,
                branch_map={dna: n2p for n2p, dna in branch_map.items()},
                limit=args.limit,
            )

    print(f"Wrote mismatch report to {output_path}")


if __name__ == "__main__":
    main()
