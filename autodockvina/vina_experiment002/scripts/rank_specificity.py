import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import os

# --- KONFIGURATION ---
# Standard-Namen (können via Argumente überschrieben werden)
DEFAULT_CSV = "/workspace/vina_experiment002/outputs/orthosteric_auto/summary_scores.csv"
OUTPUT_PREFIX = "/workspace/vina_experiment003/outputs/analysis_ranking"

def analyze_specificity(csv_file, target_receptor, reference_receptor):
    print(f"--- Lade Daten aus: {csv_file} ---")
    
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print("Fehler: CSV Datei nicht gefunden.")
        sys.exit(1)

    # 1. Datenbereinigung
    # Erstelle eindeutige ID (Library + Index)
    df['Unique_ID'] = df['Ligand_Library'] + "_" + df['Ligand_ID'].astype(str)
    # Score sicherstellen als Zahl
    df['Affinity_kcal_mol'] = pd.to_numeric(df['Affinity_kcal_mol'], errors='coerce')
    
    # Filtern auf die zwei interessanten Rezeptoren
    df_filtered = df[df['Receptor'].isin([target_receptor, reference_receptor])].copy()
    
    if df_filtered.empty:
        print(f"Fehler: Keine Daten für {target_receptor} oder {reference_receptor} gefunden.")
        print("Verfügbare Rezeptoren:", df['Receptor'].unique())
        sys.exit(1)

    # 2. Ranking berechnen (Pro Rezeptor)
    # method='min': Bei Gleichstand bekommen beide den besseren Rang
    df_filtered['Rank'] = df_filtered.groupby('Receptor')['Affinity_kcal_mol'].rank(method='min', ascending=True)

    print("Ranking berechnet. Erstelle Vergleich...")

    # 3. Pivot (Rezeptoren nebeneinander stellen)
    # Wir wollen eine Tabelle: Zeile=Ligand, Spalten=[Rank_Target, Rank_Ref, Score_Target, Score_Ref]
    pivot = df_filtered.pivot(index='Unique_ID', columns='Receptor', values=['Rank', 'Affinity_kcal_mol'])
    
    # Spalten flachklopfen (MultiIndex entfernen)
    pivot.columns = [f'{stat}_{rec}' for stat, rec in pivot.columns]
    
    # Dynamische Spaltennamen basierend auf Input
    col_rank_target = f'Rank_{target_receptor}'
    col_rank_ref = f'Rank_{reference_receptor}'
    col_score_target = f'Affinity_kcal_mol_{target_receptor}'
    col_score_ref = f'Affinity_kcal_mol_{reference_receptor}'

    # Check ob durch Pivot Zeilen verloren gingen (z.B. Ligand nur in einem Rezeptor gedockt)
    pivot = pivot.dropna()
    print(f"Vergleichbare Liganden (in beiden Rezeptoren gedockt): {len(pivot)}")

    # 4. Spezifität berechnen
    # Positive Differenz = Besserer Rank im Target (kleinere Zahl) als im Ref (große Zahl)
    # Bsp: Target Rank 1, Ref Rank 100 -> Delta = 99 (Gut!)
    pivot['Rank_Delta'] = pivot[col_rank_ref] - pivot[col_rank_target]
    
    # Sortieren: Die "spezifischsten" zuerst (Hohes Delta + Guter Target Rank)
    # Wir sortieren primär nach Delta (Spezifität), aber filtern gleich noch nach absoluter Qualität
    pivot = pivot.sort_values(by='Rank_Delta', ascending=False)

    # 5. Selektion der "Besten Spezifischen"
    # Kriterien:
    # A) Muss im Target gut binden (z.B. Top 20% des Rankings) -> Wir wollen keinen Müll, der nur zufällig spezifisch ist
    # B) Muss spezifisch sein (Rank Delta > X)
    
    top_quantile = 0.20 # Top 20%
    threshold_rank = len(pivot) * top_quantile
    
    # Filter: Guter Rank im Target UND besserer Rank im Target als im Reference
    specific_hits = pivot[
        (pivot[col_rank_target] <= threshold_rank) & 
        (pivot['Rank_Delta'] > 0)
    ]

    print(f"Gefundene spezifische Binder (Top {int(top_quantile*100)}% im Target): {len(specific_hits)}")

    # 6. Speichern
    out_csv = os.path.join(OUTPUT_PREFIX,"_specific_hits.csv")
    specific_hits.to_csv(out_csv)
    print(f"Ergebnisse gespeichert in: {out_csv}")

    # 7. Visualisierung (Scatter Plot)
    plt.figure(figsize=(10, 8))
    
    # Alle Punkte (grau)
    sns.scatterplot(data=pivot, x=col_rank_target, y=col_rank_ref, color='grey', alpha=0.3, label='All Ligands')
    
    # Spezifische Hits (rot)
    sns.scatterplot(data=specific_hits.head(50), x=col_rank_target, y=col_rank_ref, color='red', s=50, label='Top Specific Hits')

    # Diagonale Linie (y=x) -> Unspezifische Binder
    max_rank = max(pivot[col_rank_target].max(), pivot[col_rank_ref].max())
    plt.plot([0, max_rank], [0, max_rank], ls="--", c="black", label="Non-Specific (Diagonal)")

    plt.title(f"Specificity Analysis: {target_receptor} (Target) vs {reference_receptor} (Off-Target)")
    plt.xlabel(f"Rank in Target (Lower is better)")
    plt.ylabel(f"Rank in Reference (Lower is better)")
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    
    # Invertiere Achsen nicht, da 1 (bester Rank) meist links/unten intuitiv ist. 
    # Spezifische Binder sind im Bereich "Unten Links (guter Rank Target)" aber "Oben (schlechter Rank Ref)" -> Also Bereich OBEN LINKS im Plot.
    
    plt.savefig(os.path.join(OUTPUT_PREFIX,"plot.png"))
    print(f"Plot gespeichert in: {os.path.join(OUTPUT_PREFIX, '_plot.png')}")
    
    # Top 10 Ausgabe
    print("\n--- TOP 10 SPECIFIC BINDERS ---")
    print(specific_hits[[col_score_target, col_score_ref, 'Rank_Delta']].head(10))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rank and calculate specificity of ligands.")
    parser.add_argument("target", help="Exact name of the Target Receptor (as in CSV)")
    parser.add_argument("reference", help="Exact name of the Reference/Off-Target Receptor")
    parser.add_argument("--csv", default=DEFAULT_CSV, help="Path to summary_scores.csv")
    
    args = parser.parse_args()
    
    analyze_specificity(args.csv, args.target, args.reference)