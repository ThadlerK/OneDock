#!/bin/bash
set -e


# Arguments from Snakemake
COMPLEX_PDB="$1"
REF_MOL2="$2"
FINAL_DAT="$3"
LOG_FILE="$4"
STEPS_HEAT="$5"
STEPS_EQ="$6"
STEPS_PROD="$7"
STRIDE="$8"

echo "[Debug] GPU Info on this node:" >> "$LOG_FILE"
nvidia-smi >> "$LOG_FILE" 2>&1

# Derived Variables
WORK_DIR=$(dirname "$FINAL_DAT")
LOG_DIR=$(dirname "$LOG_FILE")
mkdir -p "$WORK_DIR"
mkdir -p "$LOG_DIR"

echo "[Start] Processing $COMPLEX_PDB in $WORK_DIR" >> "$LOG_FILE"

# ==============================================================================
# Step A & B: Prepare Coordinates and Parameters Together
# ==============================================================================
echo "[1/7] Separating Complex and Parametrizing..." >> "$LOG_FILE"

# 1. Extract Receptor (Everything EXCEPT LIG)
grep -v "LIG" "$COMPLEX_PDB" | grep -E "^ATOM|^HETATM" > "$WORK_DIR/receptor.pdb"

# 2. Extract Ligand (ONLY LIG)
grep "LIG" "$COMPLEX_PDB" | grep -E "^ATOM|^HETATM" > "$WORK_DIR/ligand_raw.pdb"

# 3. Convert Ligand PDB -> Mol2 using Obabel first
#    This cleans up atom names and adds connectivity (CONECT records) 
#    that Antechamber needs.
obabel -i pdb "$WORK_DIR/ligand_raw.pdb" \
       -o mol2 -O "$WORK_DIR/ligand_ready.mol2" \
       >> "$LOG_FILE" 2>&1

# 4. Run Antechamber on the COORDINATE-CORRECT Mol2
#    This generates the charges while KEEPING the atom names and coords consistent.
antechamber -i "$WORK_DIR/ligand_ready.mol2" -fi mol2 \
            -o "$WORK_DIR/lig.mol2" -fo mol2 \
            -at gaff2 -c gas -rn LIG -dr no \
            -pf y \
            >> "$LOG_FILE" 2>&1

# 5. Generate Frcomd
parmchk2 -i "$WORK_DIR/lig.mol2" -f mol2 -o "$WORK_DIR/lig.frcmod" >> "$LOG_FILE" 2>&1

# ==============================================================================
# Step C: Build Topology (tLeap)
# ==============================================================================
echo "[3/7] Running tLeap..." >> "$LOG_FILE"

cat > "$WORK_DIR/leap.in" <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams $WORK_DIR/lig.frcmod

# 1. Load the Ligand Mol2
#    CRITICAL: 'loadmol2' loads BOTH parameters AND coordinates.
#    We do NOT load a PDB for the ligand, preventing the name mismatch.
LIG = loadmol2 $WORK_DIR/lig.mol2

# 2. Load the Receptor PDB
REC = loadpdb $WORK_DIR/receptor.pdb

# 3. Combine them
COMPLEX = combine {REC LIG}

# 4. Save
#    'check' will still warn about close contacts, but it shouldn't be FATAL anymore.
check COMPLEX
saveamberparm COMPLEX $WORK_DIR/complex.prmtop $WORK_DIR/complex.inpcrd
savepdb COMPLEX $WORK_DIR/complex_check.pdb
quit
EOF

tleap -f "$WORK_DIR/leap.in" >> "$LOG_FILE" 2>&1

# D. Split Topologies (ParmEd)
echo "[4/7] Splitting Topologies..." >> "$LOG_FILE"

# OPTION 1: Explicitly remove files before running parmed (Safest/Most robust)
rm -f "$WORK_DIR/receptor.prmtop" "$WORK_DIR/ligand.prmtop"

parmed -p "$WORK_DIR/complex.prmtop" <<EOF >> "$LOG_FILE" 2>&1
strip :LIG
outparm $WORK_DIR/receptor.prmtop
quit
EOF

parmed -p "$WORK_DIR/complex.prmtop" <<EOF >> "$LOG_FILE" 2>&1
strip !:LIG
outparm $WORK_DIR/ligand.prmtop
quit
EOF

# E. Run OpenMM MD
echo "[5/7] Running OpenMM MD..." >> "$LOG_FILE"
# Set environment variables for the python script
export NSTLIM_HEAT="$STEPS_HEAT"
export NSTLIM_EQ="$STEPS_EQ"
export NSTLIM_PROD="$STEPS_PROD"

python workflow/scripts/run_openmm_md.py \
    "$WORK_DIR/complex.prmtop" \
    "$WORK_DIR/complex.inpcrd" \
    "$WORK_DIR/prod.dcd" \
    "$WORK_DIR/md.log" >> "$LOG_FILE" 2>&1

# F. Convert Trajectory
echo "[6/7] Converting Trajectory..." >> "$LOG_FILE"
cpptraj -p "$WORK_DIR/complex.prmtop" <<EOF >> "$LOG_FILE" 2>&1
trajin $WORK_DIR/prod.dcd
trajout $WORK_DIR/prod.nc netcdf
run
EOF

# G. Run MMPBSA
echo "[7/7] Running MMPBSA..." >> "$LOG_FILE"

cat > "$WORK_DIR/mmpbsa.in" <<EOF
&general
  startframe=1, endframe=99999, interval=${STRIDE},
  verbose=1, receptor_mask='!:LIG', ligand_mask=':LIG',
/
&gb
  igb=5, saltcon=0.15,
/
EOF
