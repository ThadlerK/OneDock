import sys
import os
from openmm import unit, Platform, CustomExternalForce, LangevinMiddleIntegrator
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, DCDReporter, StateDataReporter, GBn2, OBC2, HBonds, NoCutoff, PDBFile

# Inputs from Snakemake
prmtop_path = sys.argv[1]
inpcrd_path = sys.argv[2]
output_dcd  = sys.argv[3]
log_file    = sys.argv[4]

# Parameters
NSTLIM_HEAT = int(os.environ.get("NSTLIM_HEAT", "1000"))
NSTLIM_EQ   = int(os.environ.get("NSTLIM_EQ", "2000"))
NSTLIM_PROD = int(os.environ.get("NSTLIM_PROD", "5000"))

# restraints in kcal/mol/A^2 (like your sander setup)
K_HEAT = float(os.environ.get("REST_K_HEAT", "5.0"))
K_EQ   = float(os.environ.get("REST_K_EQ", "1.0"))
K_PROD = float(os.environ.get("REST_K_PROD", "0.0"))

DT = 0.001 * unit.picoseconds
TEMP = 300.0 * unit.kelvin
FRICTION = 5.0 / unit.picoseconds

print(f"[OpenMM] Loading {prmtop_path}...")
prmtop = AmberPrmtopFile(prmtop_path)
inpcrd = AmberInpcrdFile(inpcrd_path)

# OpenMM uses kJ/mol/nm^2; 1 kcal/mol/A^2 = 418.4 kJ/mol/nm^2
KCAL_PER_A2_TO_KJ_PER_NM2 = 418.4

def make_system(k_rest_kcal_a2: float, reference_positions):
    system = prmtop.createSystem(
        nonbondedMethod=NoCutoff,          # required with implicitSolvent in app-layer
        constraints=HBonds,
        implicitSolvent=OBC2
    )

    if k_rest_kcal_a2 > 0.0:
        k = k_rest_kcal_a2 * KCAL_PER_A2_TO_KJ_PER_NM2 * unit.kilojoule_per_mole / (unit.nanometer**2)
        force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addGlobalParameter("k", k)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for atom in prmtop.topology.atoms():
            if atom.name in ("N", "CA", "C", "O"):
                pos = reference_positions[atom.index]
                force.addParticle(atom.index, [pos.x, pos.y, pos.z])
        system.addForce(force)

    return system

# Platform selection
try:
    platform = Platform.getPlatformByName("CUDA")
except:
    platform = Platform.getPlatformByName("CPU")

properties = {}
if platform.getName() == "CUDA":
    properties["Precision"] = os.environ.get("OPENMM_CUDA_PRECISION", "mixed")

print("[OpenMM] Platforms available:",
      [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
print("[OpenMM] Using platform:", platform.getName(), "properties:", properties)

# Build a "base" simulation context to get reference positions and velocities
base_system = make_system(0.0, None)
integrator = LangevinMiddleIntegrator(TEMP, FRICTION, DT)
base_sim = Simulation(prmtop.topology, base_system, integrator, platform, properties)
base_sim.context.setPositions(inpcrd.positions)

# Reference positions for restraints
ref_state = base_sim.context.getState(getPositions=True)
ref_pos = ref_state.getPositions()

# --- 2. MINIMIZATION ---
print(f"[OpenMM] Using platform: {platform.getName()}")

# Minimization on restrained heat system (stable start)
print("[OpenMM] Minimizing ...")
min_system = make_system(K_HEAT, ref_pos)
min_int = LangevinMiddleIntegrator(TEMP, FRICTION, DT)
sim = Simulation(prmtop.topology, min_system, min_int, platform, properties)
sim.context.setPositions(inpcrd.positions)
sim.minimizeEnergy(maxIterations=20000)

# Check for NaN / Failure
state = sim.context.getState(getEnergy=True, getPositions=True)
pot_energy = state.getPotentialEnergy()
print(f"[DEBUG] Energy after minimization: {pot_energy}")



def run_stage(sim_in: Simulation, k_rest, nsteps, name, log_target, append_log= False, dcd_path=None, log_path=None):
    system = make_system(k_rest, ref_pos)
    integ = LangevinMiddleIntegrator(TEMP, FRICTION, DT)
    stage = Simulation(prmtop.topology, system, integ, platform, properties)

    st = sim_in.context.getState(getPositions=True, getVelocities=True)
    stage.context.setPositions(st.getPositions())
    stage.context.setVelocitiesToTemperature(50*unit.kelvin)

    if log_path:
        current_log = log_path # User provided full path (e.g. from sys.argv)

    stage.reporters.append(StateDataReporter(
        current_log, 1000,
        step=True, time=True, temperature=True, potentialEnergy=True, 
        progress=True, remainingTime=True, speed=True, totalSteps=nsteps, separator="\t",
        append=append_log
    ))
    if dcd_path:
        stage.reporters.append(DCDReporter(dcd_path, 1000))

    print(f"[OpenMM] Stage {name}: steps={nsteps}, k_rest={k_rest} kcal/mol/A^2")
    stage.step(nsteps)

    # return a sim-like object holding last state
    return stage

# Heat / Eq / Prod
sim = run_stage(sim, K_HEAT, NSTLIM_HEAT, "heat", log_path=log_file, append_log=True)
sim = run_stage(sim, K_EQ,   NSTLIM_EQ,   "eq",   log_path=log_file, append_log=True)
sim = run_stage(sim, K_PROD, NSTLIM_PROD, "prod", dcd_path=output_dcd, log_path=log_file, append_log=True)

print("[OpenMM] Done. Wrote {output_dcd}")