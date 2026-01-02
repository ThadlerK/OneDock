import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import sys
import textwrap
import argparse
import itertools
import re

# --- CONFIGURATION ---
# Base path where the receptor folders are located
# Example: /workspace/vina_experiment004/outputs/orthosteric_auto/
DEFAULT_BASE_DIR = "/workspace/vina_experiment004/outputs/orthosteric_auto"

# Folder for generated images
IMAGES_DIR_NAME = "receptor_comparison"

# Global design for plots
sns.set_theme(style="whitegrid")
# Pandas setting to prevent long SMILES from being truncated
pd.set_option('display.max_colwidth', None)

def shorten_receptor_name(name):
    """
    Shortens complex receptor names for better plotting readability.
    Example: SLC6A20_8I91_complete_occluded -> SLC6A20 (occluded)
             SLC6A19_inward_open -> SLC6A19 (inward)
    """
    if not isinstance(name, str):
        return str(name)
        
    # 1. Extract Protein Name (usually the first part before underscore)
    # We look for standard gene names like SLC6A20, ACE2, etc.
    parts = name.split('_')
    protein = parts[0]
    
    # 2. Extract Conformation/State
    # We look for keywords in the full string
    state = ""
    lower_name = name.lower()
    
    if "occluded" in lower_name:
        state = "occluded"
    elif "inward" in lower_name:
        state = "inward"
    elif "outward" in lower_name:
        state = "outward"
    elif "open" in lower_name:
        state = "open"
    elif "closed" in lower_name:
        state = "closed"
    
    # Construct short name
    if state:
        return f"{protein} ({state})"
    else:
        # Fallback: Just return the protein name if no state found
        # Or return the original if it's short enough
        if len(name) > 15:
            return protein
        return name

def load_data(base_path):
    """
    Recursively searches for CSV files in receptor subfolders
    and loads them into a single DataFrame.
    """
    print(f"Searching for results in: {base_path}")
    
    # Search pattern: Base folder -> Any subfolder (*) -> Any CSV (*.csv)
    # This covers both summary_scores.csv and RECEPTORNAME_summary.csv
    search_pattern = os.path.join(base_path, "*", "*summary.csv")
    files = glob.glob(search_pattern)
    
    if not files:
        print("WARNING: No CSV files found!")
        # Fallback attempt directly in the base folder (in case of old structure)
        search_pattern_fallback = os.path.join(base_path, "*.csv")
        files = glob.glob(search_pattern_fallback)
    
    if not files:
        print(f"ERROR: No data found under {base_path}")
        sys.exit(1)
        
    print(f"Files found: {len(files)}")
    
    df_list = []
    for f in files:
        try:
            # Load CSV
            tmp_df = pd.read_csv(f)
            df_list.append(tmp_df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
            
    if not df_list:
        print("No valid data could be loaded.")
        sys.exit(1)

    # Concatenate
    full_df = pd.concat(df_list, ignore_index=True)
    
    # Data cleaning & Types
    full_df['Affinity_kcal_mol'] = pd.to_numeric(full_df['Affinity_kcal_mol'], errors='coerce')
    
    # Create unique name if necessary
    if 'Unique_ID' not in full_df.columns:
        # Fallback: Construct ID from Library and ID
        if 'Ligand_Library' in full_df.columns and 'Ligand_ID' in full_df.columns:
            full_df['Unique_ID'] = full_df['Ligand_Library'] + "_" + full_df['Ligand_ID'].astype(str)
        elif 'Ligand' in full_df.columns:
            full_df['Unique_ID'] = full_df['Ligand']
        else:
            full_df['Unique_ID'] = "Unknown_" + full_df.index.astype(str)

    # SMILES Fallback
    if 'SMILES' not in full_df.columns:
        full_df['SMILES'] = "N/A"
        
    # --- APPLY SHORT NAMES ---
    if 'Receptor' in full_df.columns:
        # Keep original for reference if needed, but overwrite 'Receptor' for plotting simplicity
        full_df['Receptor_Original'] = full_df['Receptor']
        full_df['Receptor'] = full_df['Receptor'].apply(shorten_receptor_name)
    
    print(f"Total {len(full_df)} docking entries loaded.")
    return full_df

def plot_distribution(df, output_dir):
    """ Creates the boxplot for score distribution per library """
    print("--- Creating Distribution Plot ---")
    plt.figure(figsize=(14, 8))
    
    # Check if we have library info, otherwise use Unique_ID (shortened) or Receptor
    x_col = "Ligand_Library" if "Ligand_Library" in df.columns else "Receptor"
    
    sns.boxplot(data=df, x=x_col, y="Affinity_kcal_mol", hue="Receptor", palette="muted")

    plt.title("Score Distribution: Comparison of Ligand Libraries per Receptor")
    plt.ylabel("Affinity (kcal/mol) - Lower is better")
    plt.xlabel("Ligand Library / Group")
    plt.xticks(rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    outfile = os.path.join(output_dir, "plot_dist_by_library.png")
    plt.savefig(outfile)
    plt.close()
    print(f"Saved: {outfile}")

def plot_top_hits(df, output_dir, top_n=5):
    """ Creates barplots for the top hits of each library including SMILES """
    print("--- Creating Top Hits Plot ---")
    
    if "Ligand_Library" not in df.columns:
        print("Warning: No 'Ligand_Library' column. Skipping Top Hits per Library Plot.")
        return

    # Sort and select top N
    df_sorted = df.sort_values(by="Affinity_kcal_mol", ascending=True)
    top_hits_per_lib = df_sorted.groupby("Ligand_Library").head(top_n).copy()

    # Label Formatting (Name + SMILES)
    def format_label(row):
        smi = str(row['SMILES'])
        # Truncate after 30 characters
        if len(smi) > 30:
            smi_short = smi[:30] + "..."
        else:
            smi_short = smi
        return f"{row['Unique_ID']}\n{smi_short}"

    top_hits_per_lib['Plot_Label'] = top_hits_per_lib.apply(format_label, axis=1)

    # Plot
    g = sns.FacetGrid(top_hits_per_lib, col="Ligand_Library", col_wrap=3, height=5, aspect=1.5, sharey=False, sharex=False)
    g.map_dataframe(sns.barplot, x="Affinity_kcal_mol", y="Plot_Label", hue="Receptor", dodge=False, palette="viridis")

    g.add_legend()
    g.set_titles("Library: {col_name}")
    g.set_axis_labels("Affinity (kcal/mol)", "")
    
    # Adjust font size
    for ax in g.axes.flat:
        ax.tick_params(axis='y', labelsize=8)

    g.tight_layout()
    
    outfile = os.path.join(output_dir, "plot_top_hits_per_library.png")
    plt.savefig(outfile)
    plt.close()
    print(f"Saved: {outfile}")

def plot_heatmap(df, output_dir):
    """ Creates a heatmap (Receptor vs Ligand) """
    print("--- Creating Heatmap ---")
    
    # We use Unique_ID as column
    # Warning for too many data points
    if df['Unique_ID'].nunique() > 200:
        print("Note: Very many ligands (>200). The heatmap will be limited to the Top 50 (best scores).")
        # Filter to the top 50 global
        top_ids = df.sort_values("Affinity_kcal_mol", ascending=True)['Unique_ID'].unique()[:50]
        df_heatmap = df[df['Unique_ID'].isin(top_ids)].copy()
    else:
        df_heatmap = df.copy()

    # Pivot: Rows=Receptor, Columns=Ligand, Values=Affinity
    # For duplicates (same ligand against same receptor), take the mean
    heatmap_data = df_heatmap.pivot_table(index='Receptor', columns='Unique_ID', values='Affinity_kcal_mol', aggfunc='mean')

    # Safety Check: Is the dataframe valid for plotting?
    if heatmap_data.empty:
        print("Warning: Heatmap data is empty. Skipping heatmap generation.")
        return
        
    if heatmap_data.isnull().all().all():
        print("Warning: Heatmap data contains only NaNs. Skipping heatmap generation.")
        return

    # Dynamic size
    width = max(12, len(heatmap_data.columns) * 0.3)
    height = max(8, len(heatmap_data.index) * 0.8)
    
    plt.figure(figsize=(width, height))
    
    # Draw Heatmap (Reverse Viridis: Dark = Strong Binding/Negative Values)
    sns.heatmap(heatmap_data, annot=True, cmap='viridis_r', fmt=".1f", linewidths=.5)

    plt.title('Vina Docking Affinity (kcal/mol)', fontsize=16)
    plt.ylabel('Receptor', fontsize=12)
    plt.xlabel('Ligand', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    outfile = os.path.join(output_dir, "docking_heatmap.png")
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved: {outfile}")

def print_top_text(df):
    """ Prints the absolute top list as text """
    print("\n" + "="*60)
    print("ABSOLUTE TOP 10 HITS (with chemical structure)")
    print("="*60)
    
    df_sorted = df.sort_values(by="Affinity_kcal_mol", ascending=True)
    columns_to_show = ['Receptor', 'Unique_ID', 'Affinity_kcal_mol', 'SMILES']
    # Use only columns that exist
    valid_cols = [c for c in columns_to_show if c in df_sorted.columns]
    
    print(df_sorted[valid_cols].head(10).to_string(index=False))
    print("\n")

def perform_specificity_analysis(df, target_receptor, reference_receptor, output_dir):
    """
    Performs specificity analysis between two receptors:
    Calculates Ranks and Rank Deltas.
    Generates scatter plots (Rank vs Rank AND Score vs Score) and saves specific hits with SMILES.
    """
    # Note: target_receptor and reference_receptor strings passed here MUST match 
    # the shortened names if shorten_receptor_name was applied during load_data.
    
    print(f"--- Starting Specificity Analysis: {target_receptor} vs {reference_receptor} ---")
    
    # Filter for the two receptors of interest
    df_filtered = df[df['Receptor'].isin([target_receptor, reference_receptor])].copy()
    
    if df_filtered.empty:
        print(f"Error: No data found for {target_receptor} or {reference_receptor}.")
        # print("Available receptors:", df['Receptor'].unique()) # Avoid spamming this in loop
        return

    # --- DATA INTEGRITY CHECK ---
    # Use 'Unique_ID' as defined in load_data
    ids_target = set(df_filtered[df_filtered['Receptor'] == target_receptor]['Unique_ID'])
    ids_ref = set(df_filtered[df_filtered['Receptor'] == reference_receptor]['Unique_ID'])
    
    # Intersection
    common_ids = ids_target.intersection(ids_ref)
    
    if len(common_ids) == 0:
        print(f"ERROR: No common ligands found between {target_receptor} and {reference_receptor}!")
        return
        
    # Filter to common IDs
    df_clean = df_filtered[df_filtered['Unique_ID'].isin(common_ids)].copy()

    # 1. Calculate Ranking (Per Receptor)
    df_clean['Rank'] = df_clean.groupby('Receptor')['Affinity_kcal_mol'].rank(method='min', ascending=True)

    # 2. Pivot
    pivot = df_clean.pivot(index='Unique_ID', columns='Receptor', values=['Rank', 'Affinity_kcal_mol'])
    
    # Flatten Columns
    pivot.columns = [f'{stat}_{rec}' for stat, rec in pivot.columns]
    
    # Dynamic Column Names
    col_rank_target = f'Rank_{target_receptor}'
    col_rank_ref = f'Rank_{reference_receptor}'
    col_score_target = f'Affinity_kcal_mol_{target_receptor}'
    col_score_ref = f'Affinity_kcal_mol_{reference_receptor}'

    # Double check for NaN
    pivot = pivot.dropna()

    # 3. Calculate Specificity (Rank Delta)
    # Positive Delta = Better Rank in Target (lower number) than in Reference (higher number)
    pivot['Rank_Delta'] = pivot[col_rank_ref] - pivot[col_rank_target]
    pivot = pivot.sort_values(by='Rank_Delta', ascending=False)

    # 4. Selection (Top 20% in Target AND Positive Delta)
    top_quantile = 0.20 
    threshold_rank = len(pivot) * top_quantile
    
    specific_hits = pivot[
        (pivot[col_rank_target] <= threshold_rank) & 
        (pivot['Rank_Delta'] > 0)
    ].copy()

    # --- ADD SMILES TO RESULTS ---
    # Map Unique_ID (index of specific_hits) to SMILES from the cleaned dataframe
    if 'SMILES' in df_clean.columns:
        # We drop duplicates to get unique mapping of ID -> SMILES
        smiles_map = df_clean.drop_duplicates('Unique_ID').set_index('Unique_ID')['SMILES']
        specific_hits['SMILES'] = specific_hits.index.map(smiles_map)
    else:
        specific_hits['SMILES'] = "N/A"
    # -----------------------------

    print(f"Found {len(specific_hits)} specific binders (Top {int(top_quantile*100)}% in Target).")

    # 5. Save Results
    # Use clean filename (remove spaces/parentheses for file system safety)
    safe_target = re.sub(r'[^\w\-]', '_', target_receptor)
    safe_ref = re.sub(r'[^\w\-]', '_', reference_receptor)
    
    filename_base = f"specificity_{safe_target}_vs_{safe_ref}"
    out_csv = os.path.join(output_dir, f"{filename_base}_hits.csv")
    specific_hits.to_csv(out_csv)
    print(f"Specificity results saved to: {out_csv}")

    # 6. Scatter Plot: RANK vs RANK
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pivot, x=col_rank_target, y=col_rank_ref, color='grey', alpha=0.3, label='All Ligands')
    if not specific_hits.empty:
        sns.scatterplot(data=specific_hits.head(50), x=col_rank_target, y=col_rank_ref, color='red', s=50, label='Top Specific Hits')
    
    max_val = max(pivot[col_rank_target].max(), pivot[col_rank_ref].max())
    plt.plot([0, max_val], [0, max_val], ls="--", c="black", label="Non-Specific")

    plt.title(f"Rank Comparison: {target_receptor} vs {reference_receptor}")
    plt.xlabel(f"Rank in {target_receptor} (Lower is better)")
    plt.ylabel(f"Rank in {reference_receptor} (Lower is better)")
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plot_path_rank = os.path.join(output_dir, f"{filename_base}_plot_rank.png")
    plt.savefig(plot_path_rank)
    plt.close()
    print(f"Rank Specificity Plot saved to: {plot_path_rank}")

    # 7. Scatter Plot: SCORE vs SCORE (Affinity_kcal_mol)
    plt.figure(figsize=(10, 8))
    
    # Calculate Correlation
    correlation = pivot[col_score_target].corr(pivot[col_score_ref])
    
    sns.scatterplot(data=pivot, x=col_score_target, y=col_score_ref, color='blue', alpha=0.5, label='All Ligands')
    
    # Add diagonal line (y=x)
    min_score = min(pivot[col_score_target].min(), pivot[col_score_ref].min())
    max_score = max(pivot[col_score_target].max(), pivot[col_score_ref].max())
    plt.plot([min_score, max_score], [min_score, max_score], ls="--", c="black", label="Identical Affinity")

    plt.title(f"Affinity Comparison: {target_receptor} vs {reference_receptor}\nPearson Correlation: {correlation:.3f}")
    plt.xlabel(f"Affinity {target_receptor} (kcal/mol)")
    plt.ylabel(f"Affinity {reference_receptor} (kcal/mol)")
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plot_path_score = os.path.join(output_dir, f"{filename_base}_plot_score.png")
    plt.savefig(plot_path_score)
    plt.close()
    print(f"Score Comparison Plot saved to: {plot_path_score}")
    
    if not specific_hits.empty:
        print(specific_hits[[col_score_target, col_score_ref, 'Rank_Delta']].head(5).to_string())
    else:
        print("No specific hits found.")
    print("-" * 60)

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Docking Analysis Tool")
    parser.add_argument("--path", default=DEFAULT_BASE_DIR, help="Base directory containing receptor subfolders")
    parser.add_argument("--target", help="Name of the Target Receptor for Specificity Analysis")
    parser.add_argument("--reference", help="Name of the Reference Receptor for Specificity Analysis")
    
    args = parser.parse_args()
    
    base_dir = args.path
    images_dir = os.path.join(base_dir, IMAGES_DIR_NAME)

    # 1. Setup Directories
    if not os.path.exists(base_dir):
        print(f"Error: Base directory {base_dir} does not exist.")
        return
        
    os.makedirs(images_dir, exist_ok=True)
    
    # 2. Load Data (recursive)
    df = load_data(base_dir)
    
    # 3. Create Standard Plots
    plot_distribution(df, images_dir)
    plot_top_hits(df, images_dir)
    plot_heatmap(df, images_dir)
    
    # 4. Specificity Analysis Logic
    if args.target and args.reference:
        # Check if user input needs shortening or if they provided full names (which we might have shortened in df)
        # Strategy: Use provided args directly first. If not found, try shortening them.
        target = args.target
        ref = args.reference
        
        # Check if these exact strings exist in the DataFrame
        available_receptors = df['Receptor'].unique()
        
        if target not in available_receptors:
            short_target = shorten_receptor_name(target)
            if short_target in available_receptors:
                target = short_target
                
        if ref not in available_receptors:
            short_ref = shorten_receptor_name(ref)
            if short_ref in available_receptors:
                ref = short_ref
        
        perform_specificity_analysis(df, target, ref, images_dir)
    else:
        # No pair supplied -> All vs All
        print("\nNo specific target/reference supplied. Running ALL vs ALL specificity checks...")
        
        unique_receptors = sorted(df['Receptor'].unique())
        
        if len(unique_receptors) < 2:
            print("Not enough receptors found for specificity analysis (need at least 2).")
        else:
            pairs = list(itertools.permutations(unique_receptors, 2))
            print(f"Found {len(unique_receptors)} receptors. Processing {len(pairs)} comparison pairs.")
            
            for target, reference in pairs:
                perform_specificity_analysis(df, target, reference, images_dir)
    
    # 5. Text Output
    print_top_text(df)
    
    print(f"--- Analysis complete. All images are in: {images_dir} ---")

if __name__ == "__main__":
    main()