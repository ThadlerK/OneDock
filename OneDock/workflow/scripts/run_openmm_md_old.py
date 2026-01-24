import sys
import os
from openmm import unit, Platform, CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, DCDReporter, StateDataReporter, OBC2, HBonds, NoCutoff


# Inputs from Snakemake
prmtop_path = sys.argv[1]
inpcrd_path = sys.argv[2]
output_dcd  = sys.argv[3]
log_file    = sys.argv[4]

# Parameters
NSTLIM_HEAT = int(os.environ.get("NSTLIM_HEAT", "20000"))
NSTLIM_EQ   = int(os.environ.get("NSTLIM_EQ", "50000"))
NSTLIM_PROD = int(os.environ.get("NSTLIM_PROD", "200000"))
DT = 0.001 * unit.picoseconds
TEMP = 300.0 * unit.kelvin
FRICTION = 5.0 / unit.picoseconds

print(f"[OpenMM] Loading {prmtop_path}...")
prmtop = AmberPrmtopFile(prmtop_path)
inpcrd = AmberInpcrdFile(inpcrd_path)

ref_pos = inpcrd.positions

# Platform selection
try:
    platform = Platform.getPlatformByName("CUDA")
    props = {"Precision": "mixed"}
except:
    platform = Platform.getPlatformByName("CPU")
    props = {}

print(f"[OpenMM] Using platform: {platform.getName()}")

def make_system(k_rest_kcal_a2, reference_positions):
    # Standard Implicit Solvent Setup
    system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds, implicitSolvent=OBC2)
    
    if k_rest_kcal_a2 > 0.0:
        # 1 kcal/mol/A^2 = 418.4 kJ/mol/nm^2
        k = k_rest_kcal_a2 * 418.4 * unit.kilojoule_per_mole / (unit.nanometer**2)
        force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", k)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for atom in prmtop.topology.atoms():
            if atom.name in ("N", "CA", "C", "O"): # Backbone restraints
                pos = reference_positions[atom.index]
                force.addParticle(atom.index, [pos.x, pos.y, pos.z])
        system.addForce(force)
    return system

# 1. Minimization (Restrained)
print("[OpenMM] Minimizing...")
system_min = make_system(5.0, ref_pos)
integrator_min = LangevinMiddleIntegrator(TEMP, FRICTION, DT)
sim = Simulation(prmtop.topology, system_min, integrator_min, platform, props)
sim.context.setPositions(inpcrd.positions)

# Standard tolerance is fine if topology is correct
sim.minimizeEnergy(maxIterations=20000)

# --- DEBUG BLOCK START ---
state = sim.context.getState(getEnergy=True, getPositions=True)
pot_energy = state.getPotentialEnergy()
print(f"[DEBUG] Energy after minimization: {pot_energy}")

# 1. Check for NaN immediately
if str(pot_energy).lower() == "nan":
    print("[CRITICAL] Minimization failed! Energy is NaN.")
    sys.exit(1)

# 2. Save the minimized structure to check for distortions
from openmm.app import PDBFile
print("[DEBUG] Saving minimized.pdb...")
with open("minimized.pdb", "w") as f:
    PDBFile.writeFile(sim.topology, state.getPositions(), f)
# --- DEBUG BLOCK END ---

# Save minimized positions as reference for restraints
# state_min = sim.context.getState(getPositions=True)
# ref_pos = state_min.getPositions()

# Helper to run stages
def run_stage(prev_sim, k_rest, steps, name, write_dcd=False):
    # Re-create system with new restraint constant
    system = make_system(k_rest, ref_pos)
    integrator = LangevinMiddleIntegrator(TEMP, FRICTION, DT)
    sim_stage = Simulation(prmtop.topology, system, integrator, platform, props)
    
    # Carry over positions
    state = prev_sim.context.getState(getPositions=True)
    sim_stage.context.setPositions(state.getPositions())
    
    # CRITICAL FIX: Reset velocities to temperature to avoid instability explosion
    sim_stage.context.setVelocitiesToTemperature(50*unit.kelvin)
    
    # Reporters
    sim_stage.reporters.append(StateDataReporter(log_file, 1000, step=True, potentialEnergy=True, temperature=True, separator="\t", append=True))
    if write_dcd:
        sim_stage.reporters.append(DCDReporter(output_dcd, 1000))
        
    print(f"[OpenMM] Running {name} ({steps} steps, k={k_rest})...")
    sim_stage.step(steps)
    return sim_stage

# 2. Heat (Restrained)
sim = run_stage(sim, 5.0, NSTLIM_HEAT, "HEAT")
# 3. Equil (Less Restrained)
sim = run_stage(sim, 1.0, NSTLIM_EQ, "EQ")
# 4. Prod (Unrestrained)
sim = run_stage(sim, 0.0, NSTLIM_PROD, "PROD", write_dcd=True)

print("[OpenMM] Finished.")