#!/usr/bin/env python3
"""
Generates 1D weighted distributions of metallicities and thermodynamic states 
from Gadget snapshots. 

Plots the absolute mass fraction per bin (sum of all bins equals 1) rather than 
probability density, ensuring the y-axis is intuitively bounded between 0 and 1.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import Gadget as G

# Centralized configuration for plot styling based on the target quantity
QUANTITY_CONFIG = {
    'Z_star': {
        'xlabel': r'$\log_{10}(Z / Z_{\odot})$',
        'title_base': 'Stellar Metallicity Distribution',
        'out_base': 'stellar_metallicity'
    },
    'Z_gas': {
        'xlabel': r'$\log_{10}(Z / Z_{\odot})$',
        'title_base': 'Gas Metallicity Distribution',
        'out_base': 'gas_metallicity'
    },
    'Z_both': {
        'xlabel': r'$\log_{10}(Z / Z_{\odot})$',
        'title_base': 'Stellar vs Gas Metallicity',
        'out_base': 'metallicity_both'
    },
    'T_gas_metals': {
        'xlabel': r'$\log_{10}(T)$',
        'title_base': 'Distribution of Gas Metals by Temperature',
        'out_base': 'gas_metals_temp'
    },
    'rho_gas_metals': {
        'xlabel': r'$\log_{10}(\rho)$',
        'title_base': 'Distribution of Gas Metals by Density',
        'out_base': 'gas_metals_rho'
    }
}

def extract_data(snap, quantity):
    """
    Reads required blocks and returns (values, weights) for the targeted quantity.
    Filters out unphysical values and zero-metal particles.
    """
    try:
        if quantity == 'Z_star':
            mass = G.read_block(snap, 'iM  ', parttype=4)
            metals = G.read_block(snap, 'Zs  ', parttype=4)
            if len(mass) == 0: return None, None
            
            Zmass = metals[:, 2:].sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                Z = (Zmass / mass) / 0.02
            
            mask = Z > 0
            return np.log10(Z[mask]), mass[mask]

        elif quantity == 'Z_gas':
            mass = G.read_block(snap, 'MASS', parttype=0)
            metals = G.read_block(snap, 'Zs  ', parttype=0)
            if len(mass) == 0: return None, None
            
            Zmass = metals[:, 2:].sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                Z = (Zmass / mass) / 0.02
            
            mask = Z > 0
            return np.log10(Z[mask]), mass[mask]

        elif quantity == 'T_gas_metals':
            temp = G.read_block(snap, 'TEMP', parttype=0)
            metals = G.read_block(snap, 'Zs  ', parttype=0)
            if len(temp) == 0: return None, None
            
            Zmass = metals[:, 2:].sum(axis=1)
            mask = (temp > 0) & (Zmass > 0)
            return np.log10(temp[mask]), Zmass[mask]

        elif quantity == 'rho_gas_metals':
            rho = G.read_block(snap, 'RHO ', parttype=0)
            metals = G.read_block(snap, 'Zs  ', parttype=0)
            if len(rho) == 0: return None, None
            
            Zmass = metals[:, 2:].sum(axis=1)
            mask = (rho > 0) & (Zmass > 0)
            return np.log10(rho[mask]), Zmass[mask]

    except Exception as e:
        print(f"Error reading blocks from snapshot {snap} for {quantity}: {e}")
        return None, None

def plot_distribution(snaps, quantity='Z_star', out_name=None, bins=100, 
                      xlim=None, ylim=None, comp_style='residual', alpha=1.0,
                      title=None, labels=None):
    """Generates 1D histograms mapped explicitly by mass fraction per bin."""
    
    is_file_comparison = (len(snaps) >= 2)
    is_quant_comparison = (quantity == 'Z_both')
    is_comparison_mode = is_file_comparison or is_quant_comparison
    cfg = QUANTITY_CONFIG[quantity]

    if is_quant_comparison:
        snaps_to_read = [snaps[0], snaps[0]]
        quants_to_read = ['Z_star', 'Z_gas']
        if labels is None:
            labels = ['Stars', 'Gas']
        comp_style = 'step' 
    else:
        snaps_to_read = snaps
        quants_to_read = [quantity] * len(snaps)

    if out_name is None:
        out_name = f"{cfg['out_base']}.cmp.png" if is_comparison_mode else f"{cfg['out_base']}.png"

    if title is None:
        if is_comparison_mode:
            type_str = "residuals" if comp_style == 'residual' else "comparison"
            title = f"{' vs '.join(labels)} {type_str.capitalize()}" if labels else f"{cfg['title_base']} {type_str.capitalize()}"
        else:
            title = cfg['title_base']

    # Extract raw structural metrics
    all_vals = []
    all_weights = []
    for i in range(len(snaps_to_read)):
        v, w = extract_data(snaps_to_read[i], quants_to_read[i])
        if v is None: return
        all_vals.append(v)
        all_weights.append(w)

    if xlim is None:
        mins = [np.percentile(v, 0.1) for v in all_vals]
        maxs = [np.percentile(v, 99.9) for v in all_vals]
        xlim = [min(mins) - 0.5, max(maxs) + 0.5]

    edges = np.linspace(xlim[0], xlim[1], bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    # Calculate exact mass-fractions manually for clean bin-height scaling
    all_fractions = []
    for v, w in zip(all_vals, all_weights):
        h, _ = np.histogram(v, bins=edges, weights=w)
        total_w = np.sum(w)
        all_fractions.append(h / total_w if total_w > 0 else h)

    # --- PLOTTING ---
    plt.figure(figsize=(8, 6))

    if is_comparison_mode:
        if comp_style == 'residual':
            H1 = all_fractions[0]
            H2 = all_fractions[1]
            
            R = np.zeros(bins)
            C = np.full(bins, 'teal', dtype=object)

            both_empty = (H1 == 0) & (H2 == 0)
            H1_only    = (H1 > 0) & (H2 == 0)
            H2_only    = (H1 == 0) & (H2 > 0)
            both_valid = (H1 > 0) & (H2 > 0)

            R[both_valid] = (H1[both_valid] - H2[both_valid]) / H1[both_valid]
            R[both_empty] = np.nan
            C[both_empty] = 'white'
            R[H1_only] = 1.0
            C[H1_only] = 'gainsboro'
            R[H2_only] = -1.0
            C[H2_only] = 'gainsboro'

            plt.bar(centers, R, width=widths, color=C, edgecolor='black', alpha=0.8, align='center')
            plt.axhline(0, color='black', linewidth=1.2, linestyle='--')
            plt.ylabel('Relative Residual $(Q_1 - Q_2)/Q_1$', fontsize=14)
            
            if ylim is None:
                max_R = np.nanmax(np.abs(R[both_valid])) if np.any(both_valid) else 0.5
                y_bound = max(1.1, max_R * 1.1)
                plt.ylim(-y_bound, y_bound)
            else:
                plt.ylim(ylim[0], ylim[1])

        elif comp_style == 'step':
            color_palette = ['darkcyan', 'crimson', 'indigo', 'darkorange', 'forestgreen', 'chocolate']
            
            for i in range(len(all_fractions)):
                lbl = labels[i] if labels is not None else f'Snapshot {i+1}'
                current_alpha = 1.0 if i == 0 else alpha
                color = color_palette[i % len(color_palette)]
                
                # Using plt.stairs to plot precomputed fractions cleanly as step profiles
                plt.stairs(all_fractions[i], edges, color=color, linewidth=2.5, 
                           alpha=current_alpha, label=lbl)
            
            plt.ylabel('Mass Fraction per Bin', fontsize=14)
            plt.legend(loc='upper left', fontsize=12, frameon=True)
            if ylim: plt.ylim(ylim[0], ylim[1])

    else:
        # Standard Single Snapshot Mode using bar chart for explicit mass fractions
        plt.bar(centers, all_fractions[0], width=widths, color='teal', 
                edgecolor='black', alpha=0.8, align='center')
        plt.ylabel('Mass Fraction per Bin', fontsize=14)
        if ylim: plt.ylim(ylim[0], ylim[1])

    plt.xlim(xlim[0], xlim[1])
    plt.xlabel(cfg['xlabel'], fontsize=14)
    plt.title(title, fontsize=15, pad=15)
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')

    plt.savefig(out_name, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Successfully generated '{out_name}'")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot targeted metallicity and thermodynamic state mass fractions."
    )
    
    parser.add_argument('snap', type=str, nargs='+', 
                        help='Paths to snapshots (supports 1, 2, or N snapshots for step mode)')
    
    parser.add_argument('-q', '--quantity', type=str, default='Z_star', 
                        choices=list(QUANTITY_CONFIG.keys()),
                        help="Quantity to plot. Use 'Z_both' to plot Stars vs Gas in a single snapshot.")
    
    parser.add_argument('--comp_style', type=str, default='residual', choices=['residual', 'step'],
                        help="Comparison mode to use when multiple snapshots are provided")
    
    parser.add_argument('--alpha', type=float, default=1.0,
                        help="Transparency level [0.0 - 1.0] for overlay snapshot lines in step mode")
    
    parser.add_argument('--out', type=str, default=None, 
                        help='Name of the output plot (auto-generates if omitted)')
    parser.add_argument('--bins', type=int, default=100, help='Number of histogram bins (default: 100)')
    parser.add_argument('--xlim', type=float, nargs=2, default=None, metavar=('XMIN', 'XMAX'))
    parser.add_argument('--ylim', type=float, nargs=2, default=None, metavar=('YMIN', 'YMAX'))
    parser.add_argument('--title', type=str, default=None, help='Override automatically generated title')
    
    parser.add_argument('--labels', type=str, nargs='+', default=None, metavar='LABELS',
                        help='Labels for the snapshots (used sequentially for titles and legend strings)')
    
    args = parser.parse_args()

    # Guards
    if args.quantity == 'Z_both' and len(args.snap) > 1:
        print("Error: '-q Z_both' compares distributions within the same run and can only accept a single snapshot path.")
        sys.exit(1)

    if len(args.snap) > 2 and args.comp_style == 'residual':
        print("Error: Residual comparison mode is mathematically undefined for more than 2 snapshots.")
        sys.exit(1)        
        
    if not (0.0 <= args.alpha <= 1.0):
        print("Error: --alpha must be a float between 0.0 and 1.0")
        sys.exit(1)

    if args.labels:
        if args.quantity == 'Z_both':
            if len(args.labels) != 2:
                print("Error: --labels requires exactly 2 string arguments when plotting 'Z_both' (e.g., --labels 'Stars' 'Gas').")
                sys.exit(1)
        else:
            if len(args.snap) < 2:
                print("Error: --labels can only be applied in a multi-snapshot comparison execution.")
                sys.exit(1)
            if len(args.labels) != len(args.snap):
                print(f"Error: Number of labels ({len(args.labels)}) must match number of snapshots ({len(args.snap)}).")
                sys.exit(1)

    plot_distribution(
        snaps=args.snap, 
        quantity=args.quantity,
        out_name=args.out, 
        bins=args.bins, 
        xlim=args.xlim,
        ylim=args.ylim,
        comp_style=args.comp_style,
        alpha=args.alpha,
        title=args.title,
        labels=args.labels
    )