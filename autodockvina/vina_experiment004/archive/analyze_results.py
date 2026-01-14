import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# --- KONFIGURATION ---
file_path = "/workspace/vina_experiment004/outputs/orthosteric_occluded/summary_scores.csv"
image_output_dir = "/workspace/vina_experiment004/outputs/images"

# Sicherstellen, dass der Ausgabeordner für Bilder existiert
os.makedirs(image_output_dir, exist_ok=True)

# Pandas-Einstellung, damit lange SMILES nicht abgeschnitten werden (wichtig!)
pd.set_option('display.max_colwidth', None)

# 1. DATEN LADEN
try:
    df = pd.read_csv(file_path)
except FileNotFoundError:
    print(f"Fehler: Datei '{file_path}' nicht gefunden.")
    exit()

# Datenbereinigung
df['Affinity_kcal_mol'] = pd.to_numeric(df['Affinity_kcal_mol'], errors='coerce')

# Erstelle eindeutigen Namen: "Glutamine_0", "Glutamine_1" etc.
df['Unique_Ligand'] = df['Ligand_Library'] + "_" + df['Ligand_ID'].astype(str)

# Falls SMILES fehlen (z.B. alte CSV), füllen wir sie mit "N/A" auf, damit das Skript nicht crasht
if 'SMILES' not in df.columns:
    df['SMILES'] = "N/A"

# Globales Design
sns.set_theme(style="whitegrid")

# --- TEIL 1: DISTRIBUTION PRO LIBRARY & RECEPTOR ---
print("Erstelle Distribution Plot...")
plt.figure(figsize=(14, 8))
sns.boxplot(data=df, x="Ligand_Library", y="Affinity_kcal_mol", hue="Receptor", palette="muted")

plt.title("Score-Verteilung: Vergleich der Liganden-Libraries pro Rezeptor")
plt.ylabel("Affinity (kcal/mol) - Tiefer ist besser")
plt.xlabel("Ligand Library (Source File)")
plt.xticks(rotation=45, ha='right')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(image_output_dir, "plot_dist_by_library.png"))
plt.close()


# --- TEIL 2: TOP HITS PRO LIBRARY (PLOT) ---
print("Erstelle Top Hits Plot mit SMILES...")

top_n = 5
df_sorted = df.sort_values(by="Affinity_kcal_mol", ascending=True)
top_hits_per_lib = df_sorted.groupby("Ligand_Library").head(top_n).copy()

# TRICK: Wir bauen uns ein neues Label für die Y-Achse
# Format: "Name_ID \n SMILES (ersten 30 Zeichen...)"
def format_label(row):
    smi = str(row['SMILES'])
    # Abschneiden nach 30 Zeichen, damit der Plot nicht platzt
    if len(smi) > 30:
        smi_short = smi[:30] + "..."
    else:
        smi_short = smi
    return f"{row['Unique_Ligand']}\n{smi_short}"

top_hits_per_lib['Plot_Label'] = top_hits_per_lib.apply(format_label, axis=1)

# Plotten mit dem neuen Label
# sharey=False ist wichtig, damit jede Library ihre eigenen Labels hat
g = sns.FacetGrid(top_hits_per_lib, col="Ligand_Library", col_wrap=3, height=5, aspect=1.5, sharey=False, sharex=False)

g.map_dataframe(sns.barplot, x="Affinity_kcal_mol", y="Plot_Label", hue="Receptor", dodge=False, palette="viridis")

g.add_legend()
g.set_titles("Library: {col_name}")
g.set_axis_labels("Affinity (kcal/mol)", "") # Y-Label leer lassen, da Text selbsterklärend

# Schriftgröße der Y-Achse anpassen (damit die SMILES klein sind)
for ax in g.axes.flat:
    ax.tick_params(axis='y', labelsize=8) # Schriftgröße 8 für die Labels

g.tight_layout()
plt.savefig(os.path.join(image_output_dir, "plot_top_hits_per_library.png"))
plt.close()

print(f"Fertig! Bilder gespeichert in {image_output_dir}")

# --- TEIL 3: TEXT-OUTPUT MIT SMILES ---
print("\n" + "="*60)
print("ABSOLUTE TOP 10 HITS (mit chemischer Struktur)")
print("="*60)

# Hier wählen wir explizit die SMILES Spalte mit aus
columns_to_show = ['Receptor', 'Unique_Ligand', 'Affinity_kcal_mol', 'SMILES']

# Wir nehmen die top 10 des gesamten Datensatzes
top_10_overall = df_sorted[columns_to_show].head(10)

# Ausgabe als String ohne Index
print(top_10_overall.to_string(index=False))
print("\n")