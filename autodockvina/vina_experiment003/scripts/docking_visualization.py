import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

# --- CONFIGURATION ---
# Default path matches the bash script's structure
DEFAULT_CSV_PATH = "/workspace/vina_experiment001/data/outputs/summary_scores.csv"

def plot_results(csv_path):
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found at {csv_path}")
        print("Usage: python3 plot_docking_results.py [path_to_csv]")
        return

    print(f"Reading data from {csv_path}...")
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # Basic validation
    required_cols = {'Receptor', 'Ligand', 'Affinity_kcal_mol'}
    if not required_cols.issubset(df.columns):
        print(f"Error: CSV must contain columns: {required_cols}")
        print(f"Found: {df.columns}")
        return

    print(f"Found {len(df)} docking entries.")

    # Pivot the data for the heatmap
    # Index (Rows) = Receptor, Columns = Ligand, Values = Affinity
    try:
        heatmap_data = df.pivot(index='Receptor', columns='Ligand', values='Affinity_kcal_mol')
    except ValueError:
        print("Warning: Duplicate entries found (same receptor-ligand pair). taking the mean.")
        heatmap_data = df.pivot_table(index='Receptor', columns='Ligand', values='Affinity_kcal_mol', aggfunc='mean')

    # Setup the plot
    # Dynamic size: slightly larger if there are many ligands/receptors
    width = max(10, len(heatmap_data.columns) * 0.5)
    height = max(8, len(heatmap_data.index) * 0.5)
    plt.figure(figsize=(width, height))
    
    # Create Heatmap
    # cmap='viridis_r': Reversed Viridis. 
    # Yellow/Green = High numbers (Poor binding/Zero)
    # Purple/Dark Blue = Low numbers (Strong negative binding)
    sns.heatmap(heatmap_data, annot=True, cmap='viridis_r', fmt=".1f", linewidths=.5)

    plt.title('Vina Docking Affinity (kcal/mol)', fontsize=16)
    plt.ylabel('Receptor', fontsize=12)
    plt.xlabel('Ligand', fontsize=12)
    
    # Rotate x-axis labels if they overlap
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save
    output_dir = os.path.dirname(csv_path)
    output_file = os.path.join(output_dir, "docking_heatmap.png")
    
    plt.savefig(output_file, dpi=300)
    print(f"Heatmap saved to: {output_file}")
    
    # Attempt to show plot (might not work in headless environments)
    try:
        plt.show()
    except:
        pass

if __name__ == "__main__":
    # Allow command line argument for CSV path
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = DEFAULT_CSV_PATH
    
    plot_results(csv_file)