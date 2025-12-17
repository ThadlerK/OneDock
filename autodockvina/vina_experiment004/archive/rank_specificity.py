import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import os

# --- KONFIGURATION ---
# Pfad für Input und Output Ordner
DEFAULT_CSV = "/workspace/vina_experiment004/outputs/orthosteric_occluded/summary_scores.csv"
OUTPUT_DIR = "/workspace/vina_experiment004/outputs/analysis_ranking"

# Ordner erstellen
os.makedirs(OUTPUT_DIR, exist_ok=True)

def analyze_specificity(csv_file, target_receptor, reference_receptor):
    print(f"--- Lade Daten aus: {csv_file} ---")
    
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print("Fehler: CSV Datei nicht gefunden.")
        sys.exit(1)

    # 1. Datenbereinigung
    df['Unique_ID'] = df['Ligand_Library'] + "_" + df['Ligand_ID'].astype(str)
    df['Affinity_kcal_mol'] = pd.to_numeric(df['Affinity_kcal_mol'], errors='coerce')
    
    # Filtern auf die zwei interessanten Rezeptoren
    df_filtered = df[df['Receptor'].isin([target_receptor, reference_receptor])].copy()
    
    if df_filtered.empty:
        print(f"Fehler: Keine Daten für {target_receptor} oder {reference_receptor} gefunden.")
        print("Verfügbare Rezeptoren:", df['Receptor'].unique())
        sys.exit(1)

    # --- NEUER CHECK: ÜBERSCHNEIDUNG PRÜFEN ---
    # Wir holen alle Unique IDs pro Rezeptor
    ids_target = set(df_filtered[df_filtered['Receptor'] == target_receptor]['Unique_ID'])
    ids_ref = set(df_filtered[df_filtered['Receptor'] == reference_receptor]['Unique_ID'])
    
    # Schnittmenge (Intersection) finden
    common_ids = ids_target.intersection(ids_ref)
    
    print("-" * 40)
    print(f"CHECK DATA INTEGRITY:")
    print(f"Liganden in Target ({target_receptor}): {len(ids_target)}")
    print(f"Liganden in Reference ({reference_receptor}): {len(ids_ref)}")
    print(f"-> Liganden in BEIDEN (werden geplottet): {len(common_ids)}")
    
    if len(common_ids) == 0:
        print("FEHLER: Es gibt keine gemeinsamen Liganden zwischen den Rezeptoren!")
        print("Bitte überprüfe die Namen der Rezeptoren und ob das Docking vollständig durchlief.")
        sys.exit(1)
        
    missing_in_target = len(ids_ref) - len(common_ids)
    missing_in_ref = len(ids_target) - len(common_ids)
    if missing_in_target > 0 or missing_in_ref > 0:
        print(f"Warnung: {missing_in_target + missing_in_ref} Liganden wurden ignoriert, da sie nur in einem Datensatz vorkamen.")
    print("-" * 40)

    # Wir filtern den DataFrame jetzt schon auf die gemeinsamen IDs, um sicher zu gehen
    df_clean = df_filtered[df_filtered['Unique_ID'].isin(common_ids)].copy()

    # 2. Ranking berechnen (Pro Rezeptor)
    # method='min': Bei Gleichstand bekommen beide den besseren Rang
    df_clean['Rank'] = df_clean.groupby('Receptor')['Affinity_kcal_mol'].rank(method='min', ascending=True)

    print("Ranking berechnet. Erstelle Pivot-Tabelle...")

    # 3. Pivot (Rezeptoren nebeneinander stellen)
    pivot = df_clean.pivot(index='Unique_ID', columns='Receptor', values=['Rank', 'Affinity_kcal_mol'])
    
    # Spalten flachklopfen
    pivot.columns = [f'{stat}_{rec}' for stat, rec in pivot.columns]
    
    # Dynamische Spaltennamen
    col_rank_target = f'Rank_{target_receptor}'
    col_rank_ref = f'Rank_{reference_receptor}'
    col_score_target = f'Affinity_kcal_mol_{target_receptor}'
    col_score_ref = f'Affinity_kcal_mol_{reference_receptor}'

    # SICHERHEITS-CHECK: dropna()
    # Das entfernt Zeilen, die NaN enthalten (also wo ein Rezeptor fehlt).
    # Durch unseren Schritt oben (common_ids) sollte das hier eigentlich nichts mehr löschen,
    # aber es ist eine doppelte Absicherung.
    initial_len = len(pivot)
    pivot = pivot.dropna()
    final_len = len(pivot)
    
    if initial_len != final_len:
        print(f"Info: Weitere {initial_len - final_len} Zeilen durch dropna() entfernt (Datenlücken).")
    
    print(f"Finale Anzahl Punkte im Scatterplot: {len(pivot)}")

    # 4. Spezifität berechnen
    pivot['Rank_Delta'] = pivot[col_rank_ref] - pivot[col_rank_target]
    pivot = pivot.sort_values(by='Rank_Delta', ascending=False)

    # 5. Selektion
    top_quantile = 0.20 
    threshold_rank = len(pivot) * top_quantile
    
    specific_hits = pivot[
        (pivot[col_rank_target] <= threshold_rank) & 
        (pivot['Rank_Delta'] > 0)
    ]

    print(f"Gefundene spezifische Binder (Top {int(top_quantile*100)}% im Target): {len(specific_hits)}")

    # 6. Speichern
    out_csv = os.path.join(OUTPUT_DIR, "specific_hits.csv")
    specific_hits.to_csv(out_csv)
    print(f"Ergebnisse gespeichert in: {out_csv}")

    # 7. Visualisierung (Scatter Plot)
    plt.figure(figsize=(10, 8))
    
    # Alle Punkte (grau)
    sns.scatterplot(data=pivot, x=col_rank_target, y=col_rank_ref, color='grey', alpha=0.3, label='All Ligands')
    
    # Spezifische Hits (rot)
    if not specific_hits.empty:
        sns.scatterplot(data=specific_hits.head(50), x=col_rank_target, y=col_rank_ref, color='red', s=50, label='Top Specific Hits')

    # Diagonale Linie
    max_val = max(pivot[col_rank_target].max(), pivot[col_rank_ref].max())
    plt.plot([0, max_val], [0, max_val], ls="--", c="black", label="Non-Specific")

    plt.title(f"Specificity: {target_receptor} vs {reference_receptor}\n(Only ligands present in BOTH dockings)")
    plt.xlabel(f"Rank in Target (Lower is better)")
    plt.ylabel(f"Rank in Reference (Lower is better)")
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plot_path = os.path.join(OUTPUT_DIR, "specificity_plot.png")
    plt.savefig(plot_path)
    print(f"Plot gespeichert in: {plot_path}")
    
    # Top 10 Ausgabe
    print("\n--- TOP 10 SPECIFIC BINDERS ---")
    if not specific_hits.empty:
        # Falls SMILES im ursprünglichen DF waren, holen wir sie zurück (etwas komplexer nach Pivot)
        # Hier geben wir erstmal die Scores aus
        print(specific_hits[[col_score_target, col_score_ref, 'Rank_Delta']].head(10))
    else:
        print("Keine spezifischen Hits gefunden.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rank and calculate specificity of ligands.")
    parser.add_argument("target", help="Exact name of the Target Receptor")
    parser.add_argument("reference", help="Exact name of the Reference Receptor")
    parser.add_argument("--csv", default=DEFAULT_CSV, help="Path to summary_scores.csv")
    
    args = parser.parse_args()
    
    analyze_specificity(args.csv, args.target, args.reference)

# python3 /workspace/vina_experiment004/scripts/rank_specificity.py SLC6A20_8I91_complete_occluded SLC6A19_occluded_8I92 --csv /workspace/vina_experiment004/outputs/orthosteric_occluded/summary_scores.csv