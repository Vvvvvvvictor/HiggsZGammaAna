#!/usr/bin/env python3
"""
Script for comparing BDT scores in two ROOT files
- File 1: /eos/project-h/htozg-dy-privatemc/rzou/bdt/Output_ggF_rui_commonparam/SM1_2018_output.root: BDT_score in outtree
- File 2: /eos/home-j/jiehan/root/outputs/bdt_comparison_20250813_104019/zero_to_one_jet/ZGToLLG_2018_zero_to_one_jet_bdt_comparison.root: external_bdt_score in tree
Compare the difference between the two scores by matching the event value
"""

import os
import sys
import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import logging

# Set up logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

from pdb import set_trace

def get_args():
    """Get command line arguments"""
    parser = ArgumentParser(description='Compare BDT scores in two ROOT files')
    parser.add_argument(
        '--file1',
        type=str,
        default='/eos/project-h/htozg-dy-privatemc/rzou/bdt/Output_ggF_rui_commonparam/SM1_2018_output.root',
        help='First ROOT file path (Rui\'s result)')
    parser.add_argument('--tree1', type=str, default='outtree',
                        help='Tree name in the first file')
    parser.add_argument('--score1', type=str, default='BDT_score',
                        help='Score variable name in the first file')
    parser.add_argument('--event1', type=str, default='event',
                        help='Event number variable name in the first file')

    parser.add_argument(
        '--file2',
        type=str,
        default='/eos/user/j/jiehan/root/outputsbdt_comparison_20250815_050518/zero_to_one_jet/ZGToLLG_2018_zero_to_one_jet_bdt_comparison.root',
        help='Second ROOT file path (your result)')
    parser.add_argument('--tree2', type=str, default='tree',
                        help='Tree name in the second file')
    parser.add_argument('--score2', type=str, default='your_bdt_score',
                        help='Score variable name in the second file')
    parser.add_argument('--event2', type=str, default='event',
                        help='Event number variable name in the second file')

    parser.add_argument(
        '-o',
        '--output-dir',
        type=str,
        default='/eos/home-j/jiehan/root/outputs/bdt_score_comparison/',
        help='Output directory')
    parser.add_argument(
        '--max-events',
        type=int,
        default=None,
        help='Maximum number of events to process (for testing)')

    return parser.parse_args()


class BDTScoreComparator:
    """BDT score comparison class"""

    def __init__(self, args):
        self.args = args

        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # File and variable info
        self.file1 = args.file1
        self.tree1 = args.tree1
        self.score1 = args.score1
        self.event1 = args.event1

        self.file2 = args.file2
        self.tree2 = args.tree2
        self.score2 = args.score2
        self.event2 = args.event2

        self.output_dir = args.output_dir
        self.max_events = args.max_events

    def load_data(self):
        """Load data from two ROOT files"""
        logging.info("Start loading data...")

        # Check if files exist
        if not os.path.exists(self.file1):
            raise FileNotFoundError(f"File 1 does not exist: {self.file1}")
        if not os.path.exists(self.file2):
            raise FileNotFoundError(f"File 2 does not exist: {self.file2}")

        try:
            # Load file 1 (Rui's result)
            logging.info(f"Loading file 1: {self.file1}:{self.tree1}")
            with uproot.open(f"{self.file1}:{self.tree1}") as tree1:
                # Check if variables exist
                if self.score1 not in tree1.keys():
                    logging.error(
                        f"Variable {self.score1} does not exist in tree {self.tree1}")
                    logging.info(f"Available variables: {list(tree1.keys())}")
                    raise KeyError(f"Variable {self.score1} does not exist")

                if self.event1 not in tree1.keys():
                    logging.error(
                        f"Event number variable {self.event1} does not exist in tree {self.tree1}")
                    logging.info(f"Available variables: {list(tree1.keys())}")
                    raise KeyError(
                        f"Event number variable {self.event1} does not exist")

                # Load data
                entry_stop = self.max_events if self.max_events else None
                data1 = tree1.arrays([self.event1, self.score1, 'njet'],
                                     library="pd",
                                     entry_stop=entry_stop)
                # data1 = data1.query('njet == 0')

            logging.info(f"File 1 loaded, number of events: {len(data1)}")

            # Load file 2 (your result)
            logging.info(f"Loading file 2: {self.file2}:{self.tree2}")
            with uproot.open(f"{self.file2}:{self.tree2}") as tree2:
                # Check if variables exist
                if self.score2 not in tree2.keys():
                    logging.error(
                        f"Variable {self.score2} does not exist in tree {self.tree2}")
                    logging.info(f"Available variables: {list(tree2.keys())}")
                    raise KeyError(f"Variable {self.score2} does not exist")

                if self.event2 not in tree2.keys():
                    logging.error(
                        f"Event number variable {self.event2} does not exist in tree {self.tree2}")
                    logging.info(f"Available variables: {list(tree2.keys())}")
                    raise KeyError(
                        f"Event number variable {self.event2} does not exist")

                # Load data
                entry_stop = self.max_events if self.max_events else None
                data2 = tree2.arrays([self.event2, self.score2, "njet"],
                                     library="pd",
                                     entry_stop=entry_stop)
                # data2 = data2.query('njet == 0')

            logging.info(f"File 2 loaded, number of events: {len(data2)}")

            return data1, data2

        except Exception as e:
            logging.error(f"Error loading data: {e}")
            raise

    def match_events(self, data1, data2):
        """Match two datasets by event number"""
        logging.info("Start matching events...")

        # Rename columns for merging
        data1_renamed = data1.rename(columns={
            self.event1: 'event',
            self.score1: 'score1'
        })

        data2_renamed = data2.rename(columns={
            self.event2: 'event',
            self.score2: 'score2'
        })

        # Inner join on event number
        matched_data = pd.merge(
            data1_renamed,
            data2_renamed,
            on='event',
            how='inner')

        logging.info(f"Number of matched events: {len(matched_data)}")
        logging.info(f"Number of events in file 1: {len(data1)}")
        logging.info(f"Number of events in file 2: {len(data2)}")
        logging.info(
            f"Matching rate: {len(matched_data)/min(len(data1), len(data2))*100:.2f}%")

        if len(matched_data) == 0:
            raise ValueError(
                "No matched events! Please check if event numbers are consistent.")

        return matched_data

    def calculate_differences(self, matched_data):
        """Calculate score differences"""
        logging.info("Calculating score differences...")

        # Calculate difference: score2 - score1 (your score - Rui's score)
        matched_data['score_diff'] = matched_data['score2'] - \
            matched_data['score1']

        # Calculate relative difference
        # Avoid division by zero
        non_zero_mask = np.abs(matched_data['score1']) > 1e-10
        matched_data['score_rel_diff'] = np.nan
        matched_data.loc[non_zero_mask, 'score_rel_diff'] = (
            matched_data.loc[non_zero_mask, 'score_diff'] /
            matched_data.loc[non_zero_mask, 'score1']
        )

        # Statistics
        logging.info("=== Score Statistics ===")
        logging.info(
            f"File 1 ({self.score1}) - Mean: {matched_data['score1'].mean():.6f}, Std: {matched_data['score1'].std():.6f}")
        logging.info(
            f"File 2 ({self.score2}) - Mean: {matched_data['score2'].mean():.6f}, Std: {matched_data['score2'].std():.6f}")
        logging.info(
            f"Absolute difference - Mean: {matched_data['score_diff'].mean():.6f}, Std: {matched_data['score_diff'].std():.6f}")
        logging.info(
            f"Absolute difference - Min: {matched_data['score_diff'].min():.6f}, Max: {matched_data['score_diff'].max():.6f}")

        if not matched_data['score_rel_diff'].isna().all():
            logging.info(
                f"Relative difference - Mean: {matched_data['score_rel_diff'].mean():.6f}, Std: {matched_data['score_rel_diff'].std():.6f}")

        # Correlation
        correlation = matched_data['score1'].corr(matched_data['score2'])
        logging.info(f"Correlation: {correlation:.6f}")

        return matched_data

    def create_plots(self, matched_data):
        """Create comparison plots"""
        logging.info("Creating comparison plots...")

        # Set plot style
        plt.style.use('default')
        sns.set_palette("husl")

        # 1. Score difference histogram
        plt.figure(figsize=(12, 8))

        plt.subplot(2, 2, 1)
        plt.hist(
            matched_data['score_diff'],
            bins=50,
            alpha=0.7,
            edgecolor='black')
        plt.xlabel(f'Score Difference ({self.score2} - {self.score1})')
        plt.ylabel('Events')
        plt.title('Distribution of Score Differences')
        plt.grid(True, alpha=0.3)

        # Add statistics
        mean_diff = matched_data['score_diff'].mean()
        std_diff = matched_data['score_diff'].std()
        plt.axvline(mean_diff, color='red', linestyle='--',
                    label=f'Mean: {mean_diff:.4f}')
        plt.axvline(
            mean_diff + std_diff,
            color='orange',
            linestyle='--',
            alpha=0.7,
            label=f'±1σ: {std_diff:.4f}')
        plt.axvline(
            mean_diff -
            std_diff,
            color='orange',
            linestyle='--',
            alpha=0.7)
        plt.legend()

        # 2. Score scatter plot
        plt.subplot(2, 2, 2)
        plt.scatter(matched_data['score1'], matched_data['score2'],
                    alpha=0.6, s=1)

        # Add diagonal
        min_score = min(
            matched_data['score1'].min(),
            matched_data['score2'].min())
        max_score = max(
            matched_data['score1'].max(),
            matched_data['score2'].max())
        plt.plot([min_score, max_score], [min_score, max_score],
                 'r--', alpha=0.8, linewidth=2)

        plt.xlabel(f'{self.score1} (File 1)')
        plt.ylabel(f'{self.score2} (File 2)')
        plt.title('Score Correlation')
        plt.grid(True, alpha=0.3)

        # Add correlation
        correlation = np.corrcoef(
            matched_data['score1'],
            matched_data['score2'])[
            0,
            1]
        plt.text(0.05, 0.95, f'Correlation: {correlation:.4f}',
                 transform=plt.gca().transAxes,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # 3. Relative difference histogram (if valid data exists)
        plt.subplot(2, 2, 3)
        valid_rel_diff = matched_data['score_rel_diff'][~np.isnan(
            matched_data['score_rel_diff'])]
        if len(valid_rel_diff) > 0:
            # Remove outliers for better visualization
            percentile_99 = np.percentile(np.abs(valid_rel_diff), 99)
            filtered_rel_diff = valid_rel_diff[np.abs(
                valid_rel_diff) <= percentile_99]

            plt.hist(filtered_rel_diff, bins=50, alpha=0.7, edgecolor='black')
            plt.xlabel('Relative Difference')
            plt.ylabel('Events')
            plt.title(
                'Distribution of Relative Differences\n(99th percentile clipped)')
            plt.grid(True, alpha=0.3)

            mean_rel_diff = filtered_rel_diff.mean()
            plt.axvline(mean_rel_diff, color='red', linestyle='--',
                        label=f'Mean: {mean_rel_diff:.4f}')
            plt.legend()
        else:
            plt.text(0.5, 0.5, 'No valid relative\ndifferences to plot',
                     transform=plt.gca().transAxes, ha='center', va='center')
            plt.title('Relative Differences')

        # 4. Difference vs mean score
        plt.subplot(2, 2, 4)
        mean_scores = (matched_data['score1'] + matched_data['score2']) / 2
        plt.scatter(mean_scores, matched_data['score_diff'], alpha=0.6, s=1)
        plt.axhline(0, color='red', linestyle='--', alpha=0.8)
        plt.xlabel('Mean Score')
        plt.ylabel('Score Difference')
        plt.title('Difference vs Mean Score')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save plot
        output_path = os.path.join(self.output_dir, 'bdt_score_comparison.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        logging.info(f"Saved comparison plot: {output_path}")

        # Create detailed 2D density plot
        plt.figure(figsize=(10, 8))
        
        matched_data = matched_data.dropna(subset=['score1', 'score2', 'score_diff'])

        # 2D histogram
        plt.hist2d(matched_data['score1'], matched_data['score2'],
                   bins=50, cmap='Blues')
        plt.colorbar(label='Number of Events')

        # Add diagonal
        min_score = min(
            matched_data['score1'].min(),
            matched_data['score2'].min())
        max_score = max(
            matched_data['score1'].max(),
            matched_data['score2'].max())
        plt.plot([min_score, max_score], [min_score, max_score],
                 'r--', alpha=0.8, linewidth=2, label='y = x')

        plt.xlabel(f'{self.score1} (File 1)')
        plt.ylabel(f'{self.score2} (File 2)')
        plt.title('BDT Score Comparison - 2D Density')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Add statistics
        correlation = matched_data['score1'].corr(matched_data['score2'])
        mean_diff = matched_data['score_diff'].mean()
        std_diff = matched_data['score_diff'].std()

        info_text = f'Correlation: {correlation:.4f}\n'
        info_text += f'Mean Diff: {mean_diff:.4f}\n'
        info_text += f'Std Diff: {std_diff:.4f}\n'
        info_text += f'Events: {len(matched_data)}'

        plt.text(0.05, 0.95, info_text, transform=plt.gca().transAxes,
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                 verticalalignment='top')

        # Save 2D density plot
        output_path_2d = os.path.join(
            self.output_dir, 'bdt_score_comparison_2d.png')
        plt.tight_layout()
        plt.savefig(output_path_2d, dpi=300, bbox_inches='tight')
        plt.close()

        logging.info(f"Saved 2D density plot: {output_path_2d}")

    def save_results(self, matched_data):
        """Save results to files"""
        logging.info("Saving results...")

        # Save matched data to CSV file
        csv_path = os.path.join(self.output_dir, 'matched_scores.csv')

        # Write CSV file manually
        with open(csv_path, 'w') as f:
            f.write('event,score1,score2,score_diff,score_rel_diff\n')
            for i in range(len(matched_data['event'])):
                f.write(
                    f"{matched_data['event'][i]},{matched_data['score1'][i]},{matched_data['score2'][i]},{matched_data['score_diff'][i]},{matched_data['score_rel_diff'][i]}\n")

        logging.info(f"Saved matched data to: {csv_path}")

        # Save statistics summary to text file
        summary_path = os.path.join(self.output_dir, 'comparison_summary.txt')
        with open(summary_path, 'w') as f:
            f.write("BDT Score Comparison Summary\n")
            f.write("=" * 40 + "\n\n")

            f.write(f"File 1: {self.file1}:{self.tree1}\n")
            f.write(f"Score 1: {self.score1}\n")
            f.write(f"File 2: {self.file2}:{self.tree2}\n")
            f.write(f"Score 2: {self.score2}\n\n")

            f.write(f"Total events in file 1: {len(matched_data['event'])}\n")
            f.write(f"Total events in file 2: {len(matched_data['event'])}\n")
            f.write(f"Matched events: {len(matched_data['event'])}\n\n")

            f.write("Score Statistics:\n")
            f.write(
                f"  {self.score1} - Mean: {matched_data['score1'].mean():.6f}, Std: {matched_data['score1'].std():.6f}\n")
            f.write(
                f"  {self.score2} - Mean: {matched_data['score2'].mean():.6f}, Std: {matched_data['score2'].std():.6f}\n\n")

            f.write("Difference Statistics:\n")
            f.write(
                f"  Absolute difference - Mean: {matched_data['score_diff'].mean():.6f}, Std: {matched_data['score_diff'].std():.6f}\n")
            f.write(
                f"  Absolute difference - Min: {matched_data['score_diff'].min():.6f}, Max: {matched_data['score_diff'].max():.6f}\n")

            valid_rel_diff = matched_data['score_rel_diff'][~np.isnan(
                matched_data['score_rel_diff'])]
            if len(valid_rel_diff) > 0:
                f.write(
                    f"  Relative difference - Mean: {valid_rel_diff.mean():.6f}, Std: {valid_rel_diff.std():.6f}\n")

            correlation = np.corrcoef(
                matched_data['score1'],
                matched_data['score2'])[
                0,
                1]
            f.write(f"\nCorrelation coefficient: {correlation:.6f}\n")

        logging.info(f"Saved statistics summary to: {summary_path}")

        # Save to ROOT file
        try:
            root_path = os.path.join(self.output_dir, 'matched_scores.root')

            # Prepare data dict
            data_dict = {
                'event': matched_data['event'],
                'score1': matched_data['score1'],
                'score2': matched_data['score2'],
                'score_diff': matched_data['score_diff'],
                'score_rel_diff': matched_data['score_rel_diff']
            }

            with uproot.recreate(root_path) as f:
                f['matched_data'] = data_dict

            logging.info(f"Saved ROOT file to: {root_path}")

        except Exception as e:
            logging.warning(f"Failed to save ROOT file: {e}")

    def run(self):
        """Run comparison"""
        try:
            # Load data
            data1, data2 = self.load_data()

            # Match events
            matched_data = self.match_events(data1, data2)

            # Calculate differences
            matched_data = self.calculate_differences(matched_data)

            # Create plots
            self.create_plots(matched_data)

            # Save results
            self.save_results(matched_data)

            logging.info("Comparison finished!")

        except Exception as e:
            logging.error(f"Error during comparison: {e}")
            raise


def main():
    args = get_args()

    logging.info("Start BDT score comparison")
    logging.info(f"File 1: {args.file1}:{args.tree1}")
    logging.info(f"File 2: {args.file2}:{args.tree2}")
    logging.info(f"Output directory: {args.output_dir}")

    comparator = BDTScoreComparator(args)
    comparator.run()


if __name__ == '__main__':
    main()
