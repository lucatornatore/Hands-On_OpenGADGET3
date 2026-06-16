#!/usr/bin/env python3
"""
Creates a density-temperature phase diagram from Gadget snapshots.
Supports single snapshot profiling or dual snapshot residual calculations.
Features synchronized, optional 1D marginal mass-fraction distributions (step style).
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import matplotlib.gridspec as gridspec
import Gadget as G

def compute_density_grid(snap, Radius, xc, yc, zc, colorcode, xbins, ybins, 
                         xbins_1d=None, ybins_1d=None, plot_marginals=False):
    """Reads snapshot data and returns the 2D grid H, and optionally 1D marginals."""
    rho  = G.read_block(snap, 'RHO')
    temp = G.read_block(snap, 'TEMP') 
    mass = G.read_block(snap, 'MASS', parttype=0)
    
    if colorcode in ('metals', 'metallicity'):
        metals = G.read_block(snap, 'Zs', parttype=0)
    
    if Radius > 0:
        pos_gas = G.read_block(snap, 'POS', parttype=0)
        distance = np.linalg.norm((pos_gas - np.array([xc, yc, zc])), axis=1)
        mask = distance < Radius
        
        rho = rho[mask]
        if len(rho) == 0:
            return None
            
        temp = temp[mask]
        mass = mass[mask]
        if colorcode in ('metals', 'metallicity'):
            metals = metals[mask]

    # Calculate exactly which mass reservoir we are tracking for 1D distributions
    if colorcode in ('metals', 'metallicity'):
        weights_1d = metals[:, 2:].sum(axis=1)
    else:
        weights_1d = mass

    # 1D Marginal Mass Distributions (Normalized to mass fractions)
    H_1D_rho, H_1D_temp = None, None
    if plot_marginals and xbins_1d is not None and ybins_1d is not None:
        H_1D_rho, _  = np.histogram(rho, bins=xbins_1d, weights=weights_1d)
        H_1D_temp, _ = np.histogram(temp, bins=ybins_1d, weights=weights_1d)
        
        total_mass = np.sum(weights_1d)
        if total_mass > 0:
            H_1D_rho  = H_1D_rho / total_mass
            H_1D_temp = H_1D_temp / total_mass

    # 2D Grid Setup
    if colorcode == 'metallicity':
        Zmass = metals[:, 2:].sum(axis=1)
        H_Zmass, _, _ = np.histogram2d(rho, temp, bins=[xbins, ybins], weights=Zmass)
        H_GasMass, _, _ = np.histogram2d(rho, temp, bins=[xbins, ybins], weights=mass)
        
        with np.errstate(invalid='ignore', divide='ignore'):
            H = np.divide(H_Zmass, H_GasMass)
        H = H / 0.02
        
    else:
        if colorcode == 'mass':
            weights = mass / np.sum(mass)
        elif colorcode == 'metals':
            Zmass = metals[:, 2:].sum(axis=1)
            weights = Zmass / np.sum(Zmass)
            
        H, _, _ = np.histogram2d(rho, temp, bins=[xbins, ybins], weights=weights)
        
    return H, H_1D_rho, H_1D_temp

def compute_1d_residual(h1, h2, clim):
    """Calculates 1D residuals, capping extremes and mapping visual sentinels."""
    R = np.zeros_like(h1)
    C = np.full(len(h1), 'teal', dtype=object)
    
    both_empty = (h1 == 0) & (h2 == 0)
    h1_only    = (h1 > 0) & (h2 == 0)
    h2_only    = (h1 == 0) & (h2 > 0)
    both_valid = (h1 > 0) & (h2 > 0)

    R[both_valid] = (h1[both_valid] - h2[both_valid]) / h1[both_valid]
    
    # Cap actual numerical residuals at the colormap bounds for clean rendering
    R[both_valid] = np.clip(R[both_valid], clim[0], clim[1])

    # Assign empty bins to invisible, and missing comparisons to scale caps + gainsboro
    R[both_empty] = 0.0
    C[both_empty] = 'none' 
    R[h1_only] = clim[1]
    C[h1_only] = 'gainsboro'
    R[h2_only] = clim[0]
    C[h2_only] = 'gainsboro'
    
    return R, C

def make_phase_diagram(snaps, Radius=0.0, xc=0.0, yc=0.0, zc=0.0, 
                       colorcode='metals', clim=None, title=None,
                       xlim=None, ylim=None, plot_marginals=False):
    
    rho_first  = G.read_block(snaps[0], 'RHO')
    temp_first = G.read_block(snaps[0], 'TEMP')
    
    # 2D Bins (Fine: 300 bins)
    xbins = np.logspace(np.log10(rho_first.min()), np.log10(rho_first.max()), 301)
    ybins = np.logspace(np.log10(temp_first.min()), np.log10(temp_first.max()), 301)
    
    # 1D Bins (Coarse: 100 bins) for marginals
    xbins_1d = np.logspace(np.log10(rho_first.min()), np.log10(rho_first.max()), 101)
    ybins_1d = np.logspace(np.log10(temp_first.min()), np.log10(temp_first.max()), 101)
    
    is_residual_mode = (len(snaps) == 2)

    # Apply conditional defaults for clim before processing data
    if clim is None:
        if is_residual_mode:
            clim = [-1.0, 1.0]
        else:
            clim = [1e-5, 2.0] if colorcode == 'metallicity' else [1e-5, 5e-3]

    grid_data = compute_density_grid(snaps[0], Radius, xc, yc, zc, colorcode, 
                                     xbins, ybins, xbins_1d, ybins_1d, plot_marginals)
    if grid_data is None:
        print(f"Error: No particles found in selection for {snaps[0]}")
        return
    H1, H1_1D_rho, H1_1D_temp = grid_data

    cmap = plt.cm.get_cmap('viridis').copy()
    
    if is_residual_mode:
        grid2_data = compute_density_grid(snaps[1], Radius, xc, yc, zc, colorcode, 
                                          xbins, ybins, xbins_1d, ybins_1d, plot_marginals)
        if grid2_data is None: return
        H2, H2_1D_rho, H2_1D_temp = grid2_data
            
        plotlabel = f'Relative Residual $(Q_1 - Q_2)/Q_1$ [{colorcode}]'

        H = np.full_like(H1, np.nan)
        both_valid = (H1 > 0) & (H2 > 0)
        one_is_empty = ((H1 == 0) & (H2 > 0)) | ((H1 > 0) & (H2 == 0))

        H[both_valid] = (H1[both_valid] - H2[both_valid]) / H1[both_valid]
        H[one_is_empty] = -999.0

        cmap.set_bad(color='white', alpha=1.0)
        cmap.set_under(color='gainsboro', alpha=1.0)
        norm = cl.Normalize(vmin=clim[0], vmax=clim[1])
        
    else:
        H = H1
        label_dict = {'mass': 'mass fraction', 'metals': 'metal mass fraction', 'metallicity': 'Average Metallicity [solar]'}
        plotlabel = label_dict.get(colorcode, '.')
        cmap.set_bad(color='white', alpha=0.0) 
        norm = cl.LogNorm(vmin=clim[0], vmax=clim[1])

    # --- ADVANCED PLOTTING MATRIX ---
    if plot_marginals:
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(2, 4, width_ratios=[4, 1.2, 0.1, 0.2], height_ratios=[4, 1.2], 
                               wspace=0.05, hspace=0.05)
        
        ax_main = fig.add_subplot(gs[0, 0])
        ax_rho  = fig.add_subplot(gs[1, 0], sharex=ax_main)
        ax_temp = fig.add_subplot(gs[0, 1], sharey=ax_main)
        cax     = fig.add_subplot(gs[0, 3])
        
        # Hide intermediate ticks for clean matrix layout
        ax_main.tick_params(labelbottom=False)
        ax_temp.tick_params(labelleft=False)
    else:
        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 0.2], wspace=0.05)
        
        ax_main = fig.add_subplot(gs[0, 0])
        cax     = fig.add_subplot(gs[0, 1])
        ax_main.set_xlabel(r'$\rho$')

    # 1. 2D Phase Space Plot
    mesh = ax_main.pcolormesh(xbins, ybins, H.T, shading='auto', cmap=cmap, norm=norm)
    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_ylabel(r'$T$')
    
    # Colorbar
    fig.colorbar(mesh, cax=cax, label=plotlabel, extend='min' if is_residual_mode else 'neither')

    # 2. Marginal 1D Plot Logic
    if plot_marginals:
        widths_rho = np.diff(xbins_1d)
        heights_temp = np.diff(ybins_1d)

        if is_residual_mode:
            # Residual mode continues to use filled bars to allow for dynamic per-bin color mapping
            R_rho, C_rho = compute_1d_residual(H1_1D_rho, H2_1D_rho, clim)
            ax_rho.bar(xbins_1d[:-1], R_rho, width=widths_rho, align='edge', color=C_rho, edgecolor='none')
            ax_rho.axhline(0, color='black', lw=1.2, ls='--')
            ax_rho.set_ylim(clim[0], clim[1])
            ax_rho.set_ylabel('Residual')

            R_temp, C_temp = compute_1d_residual(H1_1D_temp, H2_1D_temp, clim)
            ax_temp.barh(ybins_1d[:-1], R_temp, height=heights_temp, align='edge', color=C_temp, edgecolor='none')
            ax_temp.axvline(0, color='black', lw=1.2, ls='--')
            ax_temp.set_xlim(clim[0], clim[1])
            ax_temp.set_xlabel('Residual')

        else:
            # Density (rho) mass fraction
            ax_rho.stairs(H1_1D_rho, xbins_1d, color='teal', linewidth=1.5, fill=False)
            ax_rho.set_yscale('log')
            ax_rho.set_ylabel('Mass Frac')

            # Temperature (T) mass fraction
            ax_temp.stairs(H1_1D_temp, ybins_1d, orientation='horizontal', color='teal', linewidth=1.5, fill=False)
            ax_temp.set_xscale('log')
            ax_temp.set_xlabel('Mass Frac')

        # Apply universal formatting for marginals
        ax_rho.set_xlabel(r'$\rho$')
        ax_rho.grid(True, linestyle='--', alpha=0.5)
        ax_temp.grid(True, linestyle='--', alpha=0.5)

    # Setting final limits forces synchronization across the shared X/Y axes
    ax_main.set_xlim(xlim[0], xlim[1])
    ax_main.set_ylim(ylim[0], ylim[1])
    
    if title: 
        fig.suptitle(title, fontsize=15, y=0.94)

    if plot_marginals:
        output_filename = f"phasespace_with_marginals.{colorcode}.cmp.png" if is_residual_mode else f"phasespace_with_marginals.{colorcode}.png"
    else:
        output_filename = f"phasespace.{colorcode}.cmp.png" if is_residual_mode else f"phasespace.{colorcode}.png"
        
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved successfully to '{output_filename}'.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate a density-temperature phase diagram with optional joint 1D distributions."
    )
    
    parser.add_argument('snap', type=str, nargs='+', 
                        help='Path to snapshot 1 (and optional snapshot 2 for residuals)')
    
    parser.add_argument('--colorcode', type=str, default='metals', choices=['metals', 'mass', 'metallicity'])
    parser.add_argument('--clim', type=float, nargs=2, default=None, metavar=('VMIN', 'VMAX'))
    parser.add_argument('--title', type=str, default=None)
    parser.add_argument('--xlim', type=float, nargs=2, default=[5e-11, 5e-3], metavar=('XMIN', 'XMAX'))
    parser.add_argument('--ylim', type=float, nargs=2, default=[5e2, 5e6], metavar=('YMIN', 'YMAX'))
    parser.add_argument('--radius', type=float, default=0.0)
    parser.add_argument('--center', type=float, nargs=3, default=[0.0, 0.0, 0.0])
    parser.add_argument('--marginals', action='store_true', help='Include 1D marginal mass distribution plots')

    args = parser.parse_args()

    if len(args.snap) > 2:
        print("Error: You can pass a maximum of 2 snapshot paths.")
        sys.exit(1)

    make_phase_diagram(
        snaps=args.snap, Radius=args.radius, xc=args.center[0], yc=args.center[1], zc=args.center[2],
        colorcode=args.colorcode, clim=args.clim, title=args.title, xlim=args.xlim, ylim=args.ylim,
        plot_marginals=args.marginals
    )
