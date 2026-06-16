#!/usr/bin/env python3
"""
plot_sfr.py

Plots star formation rate (SFR) and other thermodynamic quantities from 
cosmological simulation log files directly against Redshift (z) on the bottom axis
(with z=0 on the left), with the cosmic expansion factor (a) mapped to the top secondary axis.

Strictly non-interactive: always writes the resulting visualization directly to a PNG file.
Handles jagged/malformed log files by automatically dropping rows with column mismatches.
Includes an automatic cleanup routine to purge overwritten simulation restart histories.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
# Enforce a non-interactive backend globally to avoid loading GUI dependencies
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULT_COLUMNS = {
    1: {"name": "Expansion Factor $a$", "unit": ""},
    2: {"name": "Stochastic SFR", "unit": "code units"},
    3: {"name": "Total SFR", "unit": "code units"},
    4: {"name": "Total SFR", "unit": "$M_\\odot / \\mathrm{yr}$"},
    5: {"name": "Total Mass of Stars", "unit": "$M_\\odot$"}
}

def get_column_label(col_idx, density_mode=None, force_redshift=False):
    if force_redshift and col_idx == 1:
        return "Redshift $z$"
        
    col_info = DEFAULT_COLUMNS.get(col_idx, {"name": f"Column {col_idx}", "unit": ""})
    name = col_info["name"]
    unit = col_info["unit"]
    
    if density_mode and density_mode != "none":
        name = f"{name} Density"
        vol_unit = "$\\mathrm{Volume}^{-1}$"
        if density_mode == "physical":
            vol_unit = f"{vol_unit} (physical)"
        else:
            vol_unit = f"{vol_unit} (comoving)"
            
        if unit:
            unit = f"{unit} {vol_unit}"
        else:
            unit = vol_unit

    label = f"{name}"
    if unit:
        label += f" [{unit}]"
    return label

def a_to_z(a):
    a = np.array(a, dtype=float)
    a = np.maximum(a, 1e-10)
    return 1.0 / a - 1.0

def z_to_a(z):
    z = np.array(z, dtype=float)
    return 1.0 / (z + 1.0)

def clean_simulation_log(file_path):
    """
    Parses the log file, checks for non-monotonic behavior in the expansion factor 'a' 
    (Column 1), preserves only the final valid evolutionary timeline, renames the 
    original file to file_path.old, and writes the cleaned data to file_path.
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading file {file_path} during cleanup: {e}", file=sys.stderr)
        sys.exit(1)

    headers = [l for l in lines if l.strip().startswith('#')]
    raw_data_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]

    if not raw_data_lines:
        print(f"Error: No data found in {file_path}", file=sys.stderr)
        sys.exit(1)

    first_row_discarded = raw_data_lines[0]
    data_lines = raw_data_lines[1:]

    if not data_lines:
        return

    parsed_rows = []
    a_values = []
    
    for line in data_lines:
        parts = line.split()
        if parts:
            parsed_rows.append(line)
            try:
                a_values.append(float(parts[0]))
            except ValueError:
                # Discard non-numeric rows from the monotonic tracker
                parsed_rows.pop()

    num_entries = len(a_values)
    keep_mask = np.ones(num_entries, dtype=bool)

    max_a_from_end = -1.0
    for i in reversed(range(num_entries)):
        current_a = a_values[i]
        if max_a_from_end == -1.0:
            max_a_from_end = current_a
            continue
        
        if current_a >= max_a_from_end:
            keep_mask[i] = False
        else:
            max_a_from_end = current_a

    num_deleted = num_entries - np.sum(keep_mask)
    if num_deleted > 0:
        print(f"--> File '{file_path}' contains duplicate/overwritten states. Found {num_deleted} overlapping steps.")
        cleaned_data_lines = [first_row_discarded] + [parsed_rows[i] for i in range(num_entries) if keep_mask[i]]
        
        old_file_path = f"{file_path}.old"
        if os.path.exists(old_file_path):
            os.remove(old_file_path)
            
        os.rename(file_path, old_file_path)
        print(f"    Original file backed up as: '{old_file_path}'")
        
        with open(file_path, 'w') as f:
            if headers:
                f.writelines(headers)
            f.writelines(cleaned_data_lines)
        print(f"    Cleaned dataset successfully rewritten to: '{file_path}'")
    else:
        print(f"--> File '{file_path}' timeline is consistent and monotonic. No cleanup needed.")


def load_rectangular_data(data_lines):
    """
    Parses raw data lines and builds a rectangular matrix by skipping rows 
    that deviate from the dominant column structure of the file.
    """
    if not data_lines:
        raise ValueError("No data rows found to process.")
    
    # Split lines into whitespace-separated lists of values
    split_lines = [line.split() for line in data_lines if line.strip()]
    
    # Track the number of items in each row
    row_lengths = [len(row) for row in split_lines]
    
    # Determine the expected layout size using the statistical mode (most common width)
    # This prevents an anomaly on line 1 or the final line from biasing the dataset template.
    counts = {}
    for length in row_lengths:
        counts[length] = counts.get(length, 0) + 1
    expected_cols = max(counts, key=counts.get)
    
    cleaned_rows = []
    skipped_count = 0
    
    for idx, row in enumerate(split_lines):
        if len(row) == expected_cols:
            try:
                # Ensure all entries are numerical floats
                numerical_row = [float(val) for val in row]
                cleaned_rows.append(numerical_row)
            except ValueError:
                skipped_count += 1
        else:
            skipped_count += 1
            
    if skipped_count > 0:
        print(f"    [Warning] Filtered out {skipped_count} malformed/uneven rows to maintain layout layout uniformity ({expected_cols} columns).")
        
    return np.array(cleaned_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Plot Star Formation Rate (SFR) against Redshift (z) from cosmological logs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("files", nargs="+", help="Path to one or more simulation log files.")
    parser.add_argument("-x", "--x-col", type=int, default=1, 
                        help="1-based column index for x-axis target. If set to 1, maps to Redshift z.")
    parser.add_argument("-y1", "--y1", type=int, default=3, help="1-based column index for the left y-axis.")
    parser.add_argument("-y2", "--y2", type=int, default=None, help="1-based column index for the right y-axis.")
    parser.add_argument("--labels", nargs="+", help="Custom labels for the input files in the legend.")
    parser.add_argument("-d", "--density", action="store_true", help="Plot values as volume densities.")
    parser.add_argument("--vol-type", "--volume-type", choices=["comoving", "physical"], default="comoving",
                        dest="vol_type", help="Volume type for density calculation.")
    parser.add_argument("-V", "--volume", type=float, nargs="+", default=None, help="Cosmological simulation volume(s).")
    parser.add_argument("-s", "--point-size", type=float, default=2.0, help="Point size for the scatter plot.")
    parser.add_argument("-o", "--output", help="Output custom plot filename. If omitted, uses an automated default name.")
    parser.add_argument("--title", help="Custom title for the plot.")
    parser.add_argument("--logx", action="store_true", help="Use log scale for the x-axis (Redshift).")
    parser.add_argument("--logy1", action="store_true", help="Use log scale for the left y-axis.")
    parser.add_argument("--logy2", action="store_true", help="Use log scale for the right y-axis.")
    parser.add_argument("--xlim", type=float, nargs=2, help="X-axis limits in Redshift z (min max).")
    parser.add_argument("--ylim1", type=float, nargs=2, help="Left y-axis limits (min max).")
    parser.add_argument("--ylim2", type=float, nargs=2, help="Right y-axis limits (min max).")
    parser.add_argument("--no-expansion-top", action="store_true", help="Disable secondary top expansion factor axis.")
    parser.add_argument("--usetex", action="store_true", help="Use system LaTeX for text rendering.")
    
    args = parser.parse_args()
    
    if args.usetex:
        matplotlib.rcParams['text.usetex'] = True
    
    if args.density and args.volume is None:
        parser.error("Simulation volume (-V / --volume) is required when density mode is enabled.")
        
    if args.volume is not None:
        if len(args.volume) != 1 and len(args.volume) != len(args.files):
            parser.error(f"Number of volumes ({len(args.volume)}) must be 1 or match the number of files ({len(args.files)}).")

    if args.labels and len(args.labels) != len(args.files):
        parser.error("The number of custom labels must match the number of input files.")

    colors = ['#1f77b4', '#ff770f', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    fig, ax1 = plt.subplots(figsize=(9, 6.5))
    ax2 = ax1.twinx() if args.y2 is not None else None
        
    handles_y1 = []
    handles_y2 = []
    
    for i, file_path in enumerate(args.files):
        clean_simulation_log(file_path)

        file_label = args.labels[i] if args.labels else os.path.basename(file_path)
        color = colors[i % len(colors)]
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            data_lines = [l for l in lines if l.strip() and not l.strip().startswith('#')]
            if len(data_lines) > 0:
                data_lines = data_lines[1:] 
            else:
                raise ValueError("No data rows found in the file.")
            
            # Use the robust filtering matrix parser instead of direct np.loadtxt
            data = load_rectangular_data(data_lines)
            
        except Exception as e:
            print(f"Error reading {file_path}: {e}", file=sys.stderr)
            sys.exit(1)
            
        if len(data.shape) == 1:
            data = data.reshape(1, -1)
            
        num_cols = data.shape[1]
        
        if args.x_col > num_cols or args.x_col < 1:
            print(f"Error: X column index {args.x_col} out of bounds for {file_path}.", file=sys.stderr)
            sys.exit(1)
        if args.y1 > num_cols or args.y1 < 1:
            print(f"Error: Left y-axis column index {args.y1} out of bounds for {file_path}.", file=sys.stderr)
            sys.exit(1)
        if args.y2 is not None and (args.y2 > num_cols or args.y2 < 1):
            print(f"Error: Right y-axis column index {args.y2} out of bounds for {file_path}.", file=sys.stderr)
            sys.exit(1)
            
        raw_x = data[:, args.x_col - 1]
        y1_vals = data[:, args.y1 - 1]
        
        if args.x_col == 1:
            x_vals = a_to_z(raw_x)
        else:
            x_vals = raw_x
        
        if args.density:
            file_vol = args.volume[0] if len(args.volume) == 1 else args.volume[i]
            if args.vol_type == "comoving":
                vol = file_vol
            else:
                a_vals = data[:, 0]
                vol = file_vol * (a_vals ** 3)
            y1_plot_vals = y1_vals / vol
        else:
            y1_plot_vals = y1_vals

        lbl1 = file_label if args.y2 is None else f"{file_label} (L)"
        h1 = ax1.plot(x_vals, y1_plot_vals, marker='o', linestyle='None', 
                      color=color, markersize=args.point_size, label=lbl1)
        handles_y1.extend(h1)
        
        if args.y2 is not None:
            y2_vals = data[:, args.y2 - 1]
            y2_plot_vals = y2_vals / vol if args.density else y2_vals
                
            lbl2 = f"{file_label} (R)"
            h2 = ax2.plot(x_vals, y2_plot_vals, marker='s', linestyle='None', 
                          color=color, markersize=args.point_size, alpha=0.8, label=lbl2)
            handles_y2.extend(h2)

    if args.logx:
        ax1.set_xscale('log')
    if args.logy1:
        ax1.set_yscale('log')
    if args.logy2 and ax2 is not None:
        ax2.set_yscale('log')

    if args.xlim:
        ax1.set_xlim(args.xlim)

    if args.ylim1:
        ax1.set_ylim(args.ylim1)
    if args.ylim2 and ax2 is not None:
        ax2.set_ylim(args.ylim2)

    ax1.set_xlabel(get_column_label(args.x_col, force_redshift=True))
    ax1.set_ylabel(get_column_label(args.y1, args.vol_type if args.density else None))
    if ax2 is not None:
        ax2.set_ylabel(get_column_label(args.y2, args.vol_type if args.density else None))

    ax1.grid(True, which="both", linestyle=":", alpha=0.5)

    all_handles = handles_y1 + handles_y2
    all_labels = [h.get_label() for h in all_handles]
    
    if len(args.files) == 1 and args.y2 is not None:
        all_labels = [
            get_column_label(args.y1, args.vol_type if args.density else None),
            get_column_label(args.y2, args.vol_type if args.density else None)
        ]
        
    ax1.legend(all_handles, all_labels, loc="best", framealpha=0.9)

    if args.x_col == 1 and not args.no_expansion_top:
        try:
            secax = ax1.secondary_xaxis('top', functions=(z_to_a, a_to_z))
            secax.set_xlabel('Expansion Factor $a$')
            
            xlim_z = ax1.get_xlim()
            a_ticks = np.array([0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
            z_ticks = a_to_z(a_ticks)
            
            z_min, z_max = min(xlim_z), max(xlim_z)
            mask = (z_ticks >= z_min) & (z_ticks <= z_max)
            
            if np.any(mask):
                secax.set_ticks(a_ticks[mask])
                secax.set_xticklabels([f"{a:.2f}".rstrip('0').rstrip('.') for a in a_ticks[mask]])
        except Exception as e:
            print(f"Warning: Could not create secondary expansion axis: {e}", file=sys.stderr)

    if args.title:
        plt.title(args.title, pad=20)
    else:
        title_str = "Cosmological Simulation Star Formation History"
        if args.density:
            title_str += f" ({args.vol_type.capitalize()} Density)"
        plt.title(title_str, pad=20)

    plt.tight_layout()

    if args.output:
        output_filename = args.output
    else:
        output_filename = "sfr_vs_redshift.png" if args.x_col == 1 else "sfr_plot.png"

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved successfully to: '{output_filename}'")

if __name__ == "__main__":
    main()