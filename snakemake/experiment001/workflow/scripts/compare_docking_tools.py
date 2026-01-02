import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse
from rdkit import Chem

# --- CONFIGURATION ---
# Default paths (can be overridden by arguments)
DEFAULT_VINA_CSV = "/workspace/vina_experiment004/outputs/orthosteric_auto/receptor_comparison/specificity_SLC6A20__occluded__vs_SLC6A19__occluded__hits.csv"
DEFAULT_GLIDE_CSV = "/workspace/vina_experiment004/data/glide_results/Top20_Lib1_20_occ.csv"
OUTPUT_DIR = "/workspace/vina_experiment004/outputs/orthosteric_auto/tool_comparison"

# Glide Column Names (Adjust if needed)
GLIDE_SMILES_COL = "SMILES"
GLIDE_SCORE_COL = "r_i_docking_score"
GLIDE_TITLE_COL = "title"

# Vina Column Names
# Adjust this to match the specific score column in your specificity hits file
# Example: 'Affinity_kcal_mol_SLC6A20_8I91_complete_occluded'
# If unknown, the script will try to guess or ask user input
VINA_SCORE_COL_HINT = "Affinity_kcal_mol_SLC6A20 (occluded)" 

def canonicalize_smiles(smi):
    """ Converts SMILES to a unique canonical form """
    if not isinstance(smi, str) or pd.isna(smi):
        return None
    try:
        # Glide SMILES often contain extra info like title separated by : or space
        # We take the first part if it looks like a compound string
        clean_smi = smi.split()[0].split(':')[0]
        
        mol = Chem.MolFromSmiles(clean_smi)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except:
        return None
    return None

def main():
    parser = argparse.ArgumentParser(description="Visualize comparison between Vina and Glide results.")
    parser.add_argument("--vina", default=DEFAULT_VINA_CSV, help="Path to Vina CSV (specificity hits)")
    parser.add_argument("--glide", default=DEFAULT_GLIDE_CSV, help="Path to Glide CSV")
    parser.add_argument("--out", default=OUTPUT_DIR, help="Output directory")
    
    args = parser.parse_args()
    
    os.makedirs(args.out, exist_ok=True)
    
    print("--- Loading Data ---")
    
    # 1. Load Vina
    try:
        df_vina = pd.read_csv(args.vina)
        print(f"Vina entries: {len(df_vina)}")
    except Exception as e:
        print(f"Error loading Vina CSV: {e}")
        return

    # 2. Load Glide
    try:
        # Try semicolon separator first (common for Maestro exports)
        try:
            df_glide = pd.read_csv(args.glide, sep=';')
            if GLIDE_SMILES_COL not in df_glide.columns:
                 df_glide = pd.read_csv(args.glide, sep=',')
        except:
            df_glide = pd.read_csv(args.glide)
            
        print(f"Glide entries: {len(df_glide)}")
        
        if GLIDE_SMILES_COL not in df_glide.columns:
            print(f"Error: Glide CSV missing column '{GLIDE_SMILES_COL}'. Columns found: {list(df_glide.columns)}")
            return
            
    except Exception as e:
        print(f"Error loading Glide CSV: {e}")
        return

    # 3. Canonicalize
    print("\n--- Canonicalizing Structures ---")
    if 'SMILES' not in df_vina.columns:
        print("Error: Vina CSV must contain a 'SMILES' column.")
        return

    df_vina['Canon_SMILES'] = df_vina['SMILES'].apply(canonicalize_smiles)
    df_glide['Canon_SMILES'] = df_glide[GLIDE_SMILES_COL].apply(canonicalize_smiles)

    # Drop Invalid
    df_vina = df_vina.dropna(subset=['Canon_SMILES'])
    df_glide = df_glide.dropna(subset=['Canon_SMILES'])

    # 4. Merge
    print("\n--- Merging Results ---")
    merged = pd.merge(
        df_vina, 
        df_glide, 
        on='Canon_SMILES', 
        suffixes=('_Vina', '_Glide')
    )

    print(f"Common ligands found: {len(merged)}")

    if len(merged) == 0:
        print("No matches found based on canonical SMILES.")
        return

    # 5. Identify Vina Score Column
    vina_score_col = None
    possible_cols = [c for c in merged.columns if VINA_SCORE_COL_HINT in c and "_Vina" not in c]
    
    if len(possible_cols) == 1:
        vina_score_col = possible_cols[0]
    elif len(possible_cols) > 1:
        vina_score_col = possible_cols[0]
        print(f"Warning: Multiple score columns found. Using '{vina_score_col}'.")
    else:
        print("Error: Could not identify Vina score column.")
        return

    print(f"Comparing Vina Column: '{vina_score_col}' vs Glide Column: '{GLIDE_SCORE_COL}'")

    # --- 6A. Statistics & Plotting (SCORE) ---
    correlation = merged[vina_score_col].corr(merged[GLIDE_SCORE_COL])
    print(f"Pearson Correlation (Score): {correlation:.3f}")

    plt.figure(figsize=(10, 8))
    sns.set_theme(style="whitegrid")
    
    sns.scatterplot(data=merged, x=vina_score_col, y=GLIDE_SCORE_COL, alpha=0.7, color='blue', s=80)
    sns.regplot(data=merged, x=vina_score_col, y=GLIDE_SCORE_COL, scatter=False, color="red")

    plt.title(f"Docking Tool Comparison: Vina vs Glide (Scores)\nCorrelation: {correlation:.3f}")
    plt.xlabel(f"Vina Score ({vina_score_col})")
    plt.ylabel(f"Glide Score ({GLIDE_SCORE_COL})")
    plt.grid(True, linestyle=":", alpha=0.6)

    plot_path_score = os.path.join(args.out, "comparison_plot_score.png")
    plt.savefig(plot_path_score)
    plt.close()
    print(f"Score Plot saved to: {plot_path_score}")

    # --- 6B. Statistics & Plotting (RANK) ---
    merged['Rank_Vina'] = merged[vina_score_col].rank(method='min', ascending=True)
    merged['Rank_Glide'] = merged[GLIDE_SCORE_COL].rank(method='min', ascending=True)

    spearman_corr = merged[vina_score_col].corr(merged[GLIDE_SCORE_COL], method='spearman')
    print(f"Spearman Correlation (Rank): {spearman_corr:.3f}")

    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=merged, x='Rank_Vina', y='Rank_Glide', alpha=0.7, color='green', s=80)
    
    max_rank = max(merged['Rank_Vina'].max(), merged['Rank_Glide'].max())
    plt.plot([0, max_rank], [0, max_rank], ls="--", c="black", label="Identical Rank")

    plt.title(f"Docking Tool Comparison: Vina vs Glide (Ranks)\nSpearman Correlation: {spearman_corr:.3f}")
    plt.xlabel(f"Vina Rank")
    plt.ylabel(f"Glide Rank")
    plt.legend()
    plt.grid(True, linestyle=":", alpha=0.6)

    plot_path_rank = os.path.join(args.out, "comparison_plot_rank.png")
    plt.savefig(plot_path_rank)
    plt.close()
    print(f"Rank Plot saved to: {plot_path_rank}")

    # --- 6C. Statistics & Plotting (DISTRIBUTION) ---
    print("\n--- Generating Distribution Plot ---")
    plt.figure(figsize=(10, 6))
    
    # Plotting Vina Distribution
    sns.histplot(
        data=merged, 
        x=vina_score_col, 
        color="skyblue", 
        label="Vina", 
        kde=True, 
        element="step", 
        alpha=0.5
    )
    
    # Plotting Glide Distribution
    sns.histplot(
        data=merged, 
        x=GLIDE_SCORE_COL, 
        color="orange", 
        label="Glide", 
        kde=True, 
        element="step", 
        alpha=0.5
    )

    plt.title("Distribution of Docking Scores: Vina vs Glide")
    plt.xlabel("Docking Score")
    plt.ylabel("Frequency / Count")
    plt.legend()
    plt.grid(True, linestyle=":", alpha=0.6)

    plot_path_dist = os.path.join(args.out, "comparison_plot_distribution.png")
    plt.savefig(plot_path_dist)
    plt.close()
    print(f"Distribution Plot saved to: {plot_path_dist}")

    # 7. Save Top Hits Comparison
    top_merged = merged.sort_values(by=vina_score_col, ascending=True).head(20)
    
    print("\n--- Top 10 Vina Hits vs Glide ---")
    cols_display = ['Unique_Ligand' if 'Unique_Ligand' in top_merged.columns else 'Ligand_ID', 'Rank_Vina', vina_score_col, 'Rank_Glide', GLIDE_SCORE_COL]
    cols_display = [c for c in cols_display if c in top_merged.columns]
    
    print(top_merged[cols_display].to_string(index=False))

    csv_path = os.path.join(args.out, "comparison_merged_data.csv")
    merged.to_csv(csv_path, index=False)
    print(f"\nMerged dataset saved to: {csv_path}")

if __name__ == "__main__":
    main()