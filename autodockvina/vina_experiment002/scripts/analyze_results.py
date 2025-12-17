import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. DATEN LADEN
file_path = "/workspace/vina_experiment002/outputs/orthosteric_auto/summary_scores.csv" 
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Fehler: Datei '{file_path}' nicht gefunden.")
    exit()

# Datenbereinigung
df['Affinity_kcal_mol'] = pd.to_numeric(df['Affinity_kcal_mol'], errors='coerce')
# Erstelle eindeutigen Namen: "Glutamine_0", "Glutamine_1" etc.
df['Unique_Ligand'] = df['Ligand_Library'] + "_" + df['Ligand_ID'].astype(str)

# Globales Design
sns.set_theme(style="whitegrid")

# --- TEIL 1: DISTRIBUTION PRO LIBRARY & RECEPTOR ---
# Frage: Welche Library liefert generell die besten Scores? 
# Und gibt es Unterschiede zwischen den Rezeptoren?

plt.figure(figsize=(14, 8))
# x = Die Library (SML Datei), y = Score, hue = Rezeptor (Farbe)
sns.boxplot(data=df, x="Ligand_Library", y="Affinity_kcal_mol", hue="Receptor", palette="muted")

plt.title("Score-Verteilung: Vergleich der Liganden-Libraries pro Rezeptor")
plt.ylabel("Affinity (kcal/mol) - Tiefer ist besser")
plt.xlabel("Ligand Library (Source File)")
plt.xticks(rotation=45, ha='right') # Labels schräg stellen für Lesbarkeit
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') # Legende nach außen
plt.tight_layout()
plt.savefig("/workspace/vina_experiment002/outputs/images/plot_dist_by_library.png")
plt.close()


# --- TEIL 2: TOP HITS PRO LIBRARY ---
# Wir wollen die Top 5 besten Liganden aus JEDER Library sehen.

# Sortieren (beste zuerst) und dann gruppieren
top_n = 5
df_sorted = df.sort_values(by="Affinity_kcal_mol", ascending=True)
top_hits_per_lib = df_sorted.groupby("Ligand_Library").head(top_n)

# FacetGrid: Erstellt für jede Library einen eigenen kleinen Plot
g = sns.FacetGrid(top_hits_per_lib, col="Ligand_Library", col_wrap=3, height=4, sharey=False, sharex=False)
g.map_dataframe(sns.barplot, x="Affinity_kcal_mol", y="Ligand_ID", hue="Receptor", dodge=False, palette="viridis")

g.add_legend()
g.set_titles("Library: {col_name}")
g.set_axis_labels("Affinity (kcal/mol)", "Top Candidates")
g.tight_layout()
plt.savefig("/workspace/vina_experiment002/outputs/images/plot_top_hits_per_library.png")
plt.close()

print("Fertig! Bilder gespeichert:")
print("1. plot_dist_by_library.png (Welche Library ist generell gut?)")
print("2. plot_top_hits_per_library.png (Die besten 5 aus jeder Library)")

# Kurzer Text-Output der allerbesten Hits
print("\n--- Absolute Top 10 Hits (Library übergreifend) ---")
print(df_sorted[['Receptor', 'Ligand_ID', 'Affinity_kcal_mol']].head(10).to_string(index=False))