import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

# ==========================================
# CONFIGURATION
# ==========================================
# Define output directory for plots
output_dir = "/workspace/mid_term_presentation_plots/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# File Paths
# 1. Glide Data Path
glide_file_path = "/workspace/mid_term_presentation_plots/Chembl_19_20_occ_data_sorted.csv"

# 2. Vina Data Paths
# Note: To compare two receptors, we typically need to load two result files.
# I have filled in the first one based on your request. 
vina_file_path_1 = "/workspace/vina_experiment004/outputs/orthosteric_auto/SLC6A20_8I91_complete_occluded/SLC6A20_8I91_complete_occluded_summary.csv"
vina_file_path_2 = "/workspace/vina_experiment004/outputs/orthosteric_auto/SLC6A19_occluded_8I92/SLC6A19_occluded_8I92_summary.csv" 

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def get_top_specific_hits(df, col_target, col_ref, id_col, top_n=3):
    """
    Identifies top N ligands specific to the target receptor.
    Specificity is determined by the 'Delta Rank':
    Delta Rank = Rank(Reference) - Rank(Target)
    
    Logic:
    - Rank 1 = Best Affinity (Lowest Score).
    - We want Low Rank for Target (e.g., 1) and High Rank for Reference (e.g., 100).
    - Maximize (Rank_Ref - Rank_Target).
    """
    # Work on a copy and drop missing values to ensure fair ranking
    sub = df.dropna(subset=[col_target, col_ref]).copy()
    
    # Rank Ascending (Lowest score = Rank 1, which is "High Affinity")
    sub['rank_target'] = sub[col_target].rank(ascending=True)
    sub['rank_ref'] = sub[col_ref].rank(ascending=True)
    
    # Calculate Delta Rank
    sub['delta_rank'] = sub['rank_ref'] - sub['rank_target']
    
    # Sort by Delta Rank (Descending) to find most specific to Target
    sub_sorted = sub.sort_values('delta_rank', ascending=False)
    
    return sub_sorted.head(top_n)

def get_top_specific_hits_2(df, col_target, col_ref, id_col, top_n=3):
    """
    1. Filters for high affinity binders on the target (Top 20% quantile).
    2. Sorts the survivors by the biggest gap in Score (Thermodynamic Selectivity).
    """
    # 1. Clean data
    sub = df.dropna(subset=[col_target, col_ref]).copy()
    
    # 2. FILTER: Keep only 'Good Binders' for the target
    # Calculate the score threshold (e.g., top 20th percentile)
    # Lower score is better, so we want the quantile 0.2
    threshold = sub[col_target].quantile(0.20) 
    
    # Keep only rows where Target Score is better (lower) than threshold
    high_affinity_sub = sub[sub[col_target] <= threshold].copy()
    
    # Fallback: If filter is too strict and empties the list, use original
    if high_affinity_sub.empty:
        high_affinity_sub = sub.copy()

    # 3. SCORE DIFFERENCE (Delta Delta G)
    # Positive value means Ref score is worse (higher) than Target (lower)
    # Example: Ref (-5.0) - Target (-10.0) = +5.0 (Good Selectivity)
    high_affinity_sub['delta_score'] = high_affinity_sub[col_ref] - high_affinity_sub[col_target]
    
    # 4. Sort by Delta Score descending
    sub_sorted = high_affinity_sub.sort_values('delta_score', ascending=False)
    
    return sub_sorted.head(top_n)

def annotate_top_hits(ax, top_hits, x_col, y_col, id_col, x_offset=0.2, y_offset=0.2):
    """
    Annotates the points on the plot.
    """
    # Simple collision avoidance for labels by cycling offsets
    offsets = [(20, -5), (20, -12), (-30, -50)] 

    for i, (idx, row) in enumerate(top_hits.iterrows()):
        x_val = row[x_col]
        y_val = row[y_col]
        label = f"({x_val:.1f} / {y_val:.1f})"
        
        # Highlight the point
        ax.scatter(x_val, y_val, color='red', s=60, edgecolor='black', zorder=10)
        
        # Add Text
        xytext = offsets[i % len(offsets)]
        ax.annotate(label, (x_val, y_val), 
                    xytext=xytext, textcoords='offset points',
                    fontsize=9, fontweight='bold', color='darkred',
                    arrowprops=dict(arrowstyle="->", color='gray', connectionstyle="arc3,rad=0.2"))

# ==========================================
# PLOTTING FUNCTIONS
# ==========================================
def plot_glide_comparison(df, filename):
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # 1. Split data based on source
    # 'beide' means it docked to both. Anything else means it failed on one (score=5).
    mask_both = df['source'] == 'beide'
    df_both = df[mask_both]
    df_single = df[~mask_both]
    
    # 2. Plot "Single" receptor hits (Grey)
    if not df_single.empty:
        sns.scatterplot(
            data=df_single, 
            x='SLC6A19_score', 
            y='SLC6A20_score', 
            color='grey', 
            alpha=0.6, 
            s=60, 
            label='Docked only one Receptor'
        )
        
    # 3. Plot "Both" receptor hits (Blue)
    # We calculate correlation only on the valid 'beide' data
    correlation = np.nan
    if not df_both.empty:
        correlation = df_both['SLC6A19_score'].corr(df_both['SLC6A20_score'])
        
        sns.scatterplot(
            data=df_both, 
            x='SLC6A19_score', 
            y='SLC6A20_score', 
            color='blue', 
            alpha=0.7, 
            s=60, 
            label='Docked to Both'
        )
        # --- RANKING & HIGHLIGHTING ---
        # Target: SLC6A20 (y_axis), Reference: SLC6A19 (x_axis)
        top_3 = get_top_specific_hits_2(df_both, col_target='SLC6A20_score', col_ref='SLC6A19_score', id_col='SMILES')
        annotate_top_hits(ax, top_3, 'SLC6A19_score', 'SLC6A20_score', 'SMILES')

    # 4. Diagonal & Decor
    # if not df.empty:
    #     all_vals = pd.concat([df['SLC6A19_score'], df['SLC6A20_score']])
    #     min_v, max_v = all_vals.min(), all_vals.max()
    #     plt.plot([min_v, max_v], [min_v, max_v], 'k--', alpha=0.5, label="x=y")

    # Labels and Title
    title_text = "Glide Docking Score Comparison"
    if not np.isnan(correlation):
        title_text += f"\nPearson Correlation: {correlation:.3f}"
        
    plt.title(title_text, fontsize=14)
    plt.xlabel("Affinity SLC6A19 (kcal/mol)", fontsize=12)
    plt.ylabel("Affinity SLC6A20 (kcal/mol)", fontsize=12)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Glide plot saved to: {save_path}")

def plot_vina_comparison(df, x_col, y_col, filename):
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # Drop NaNs
    df_clean = df.dropna(subset=[x_col, y_col]).reset_index()
    
    if df_clean.empty:
        print("No valid overlapping Vina data found.")
        return

    # Correlation
    correlation = df_clean[x_col].corr(df_clean[y_col])
    
    # Scatter Plot
    sns.scatterplot(
        data=df_clean, 
        x=x_col, 
        y=y_col, 
        color='blue', 
        alpha=0.6, 
        s=60,
        label='Similar Ligands'
    )

    # --- RANKING & HIGHLIGHTING ---
    # Target: SLC6A20. Need to find which column corresponds to SLC6A20.
    # Assuming y_col is SLC6A20 based on previous context, but let's check strings.
    target_col = y_col if "SLC6A20" in y_col else x_col
    ref_col = x_col if target_col == y_col else y_col

    print(f"Ranking Vina: Target (SLC6A20) = {target_col}, Reference = {ref_col}")
    
    top_3 = get_top_specific_hits_2(df_clean, col_target=target_col, col_ref=ref_col, id_col='Ligand_ID')
    annotate_top_hits(ax, top_3, x_col, y_col, 'Ligand_ID')
    
    # Diagonal
    # min_v = min(df_clean[x_col].min(), df_clean[y_col].min())
    # max_v = max(df_clean[x_col].max(), df_clean[y_col].max())
    # plt.plot([min_v, max_v], [min_v, max_v], 'k--', alpha=0.5, label="x=y")
    
    # Labels
    x_col_parts = x_col.split('_')
    x_col_name = x_col_parts[0]
    y_col_parts = y_col.split('_')
    y_col_name = y_col_parts[0]

    plt.title(f"Vina Affinity Comparison\nPearson Correlation: {correlation:.3f}", fontsize=14)
    plt.xlabel(f"Affinity {x_col_name} (kcal/mol)", fontsize=12)
    plt.ylabel(f"Affinity {y_col_name} (kcal/mol)", fontsize=12)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Vina plot saved to: {save_path}")

# ==========================================
# MAIN EXECUTION
# ==========================================

# --- 1. GLIDE ---
if os.path.exists(glide_file_path):
    print("Processing Glide Data...")
    try:
        # Load Glide data (separator is ';')
        glide_data = pd.read_csv(glide_file_path, sep=';')
        plot_glide_comparison(glide_data, "glide_comparison_plot.png")
    except Exception as e:
        print(f"Error processing Glide data: {e}")
else:
    print(f"Glide file not found: {glide_file_path}")

# --- 2. VINA ---
print("Processing Vina Data...")
dfs = []
# Load first file
if os.path.exists(vina_file_path_1):
    dfs.append(pd.read_csv(vina_file_path_1))
else:
    print(f"Vina file 1 not found: {vina_file_path_1}")

# Load second file
if os.path.exists(vina_file_path_2):
    dfs.append(pd.read_csv(vina_file_path_2))
else:
    print(f"Vina file 2 not found: {vina_file_path_2}")

if len(dfs) == 2:
    try:
        # Merge and Pivot
        vina_combined = pd.concat(dfs, ignore_index=True)
        # Pivot: Index=Ligand_ID, Columns=Receptor, Values=Affinity
        vina_pivot = vina_combined.pivot_table(index='Ligand_ID', columns='Receptor', values='Affinity_kcal_mol')
        
        # Get column names (Receptors)
        receptors = vina_pivot.columns.tolist()
        if len(receptors) >= 2:
            col_19 = next((c for c in receptors if "SLC6A19" in c), receptors[0])
            col_20 = next((c for c in receptors if "SLC6A20" in c), receptors[1])

            plot_vina_comparison(vina_pivot, col_19, col_20, "vina_comparison_plot.png")
        else:
            print(f"Error: Found less than 2 receptors in Vina files: {receptors}")
    except Exception as e:
        print(f"Error processing Vina data: {e}")
else:
    print("Skipping Vina plot: Need both Vina files to exist.")