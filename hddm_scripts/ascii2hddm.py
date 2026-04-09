#!/usr/bin/env python3
"""
ascii2hddm.py - ASCII to HDDM converter with flexible final states

Automatically parses final state string to determine particle types and column structure.
Supports any combination of particles in any order.

Usage: ascii2hddm.py <final_state> <target> <input> <output> [options]

Examples:
  ascii2hddm.py epem p input.txt output.hddm
  ascii2hddm.py em p single_electron.txt output.hddm
  ascii2hddm.py pi0pippim p three_pion.txt output.hddm
  ascii2hddm.py ememem p triple_electron.txt output.hddm
  ascii2hddm.py mupmumpippim pb208 mixed.txt output.hddm
  ascii2hddm.py epem p input.txt output.hddm --run 30274
  ascii2hddm.py pi0 p input.txt output.hddm --run 30274 --vertex '0:0:50:80'

The script automatically:
  - Parses the final state string into individual particles
  - Determines expected column count (beam + 4*N_particles + 4*recoil)
  - Creates correct number of products in HDDM
  - Assigns particles in the order specified
"""

import hddm_s
import sys
import argparse
import random
import numpy as np
import math
import re

def fix_exponent_format(line):
    """Insert 'E' between digits and +/− for Fortran-style exponential notation"""
    fixed = ""
    prev = ""
    for char in line:
        if (char in "+-" and prev.isdigit()):
            fixed += "E"
        fixed += char
        prev = char
    return fixed


def parse_final_state(final_state_string):
    """
    Parse a final state string into individual particle names.
    
    Examples:
        "epem" -> ["em", "ep"]
        "pi0" -> ["pi0"]
        "pi0pippim" -> ["pi0", "pip", "pim"]
        "ememem" -> ["em", "em", "em"]
        "mupmum" -> ["mup", "mum"]
    
    Returns:
        List of particle names in order they appear
    """
    # Define all known particle tokens (order matters - check longer strings first!)
    particle_tokens = [
        'pi0',   # Must come before 'pip', 'pim'
        'pip', 'pim',
        'mup', 'mum',
        'ep', 'em',
        'kp', 'km',   # For future kaon support
        'eta',        # For future eta support
    ]
    
    particles = []
    remaining = final_state_string.lower()
    
    while remaining:
        matched = False
        for token in particle_tokens:
            if remaining.startswith(token):
                particles.append(token)
                remaining = remaining[len(token):]
                matched = True
                break
        
        if not matched:
            # Could not parse - provide helpful error
            print(f"Error: Could not parse '{remaining}' in final state '{final_state_string}'")
            print(f"Known particles: {', '.join(particle_tokens)}")
            sys.exit(1)
    
    return particles


def get_particle_info(particle_name):
    """
    Get particle properties from name.
    
    Returns:
        dict with keys: name, geant_pid, pdg_type, mass, charge
    """
    particle_db = {
        'photon':  {'geant_pid': 1,  'pdg_type': 22,         'mass': 0,            'charge': 0},
        'em':      {'geant_pid': 3,  'pdg_type': 11,         'mass': 0.000510999,  'charge': -1},
        'ep':      {'geant_pid': 2,  'pdg_type': -11,        'mass': 0.000510999,  'charge': 1},
        'mum':     {'geant_pid': 6,  'pdg_type': 13,         'mass': 0.105658389,  'charge': -1},
        'mup':     {'geant_pid': 5,  'pdg_type': -13,        'mass': 0.105658389,  'charge': 1},
        'pim':     {'geant_pid': 9,  'pdg_type': -211,       'mass': 0.13956995,   'charge': -1},
        'pip':     {'geant_pid': 8,  'pdg_type': 211,        'mass': 0.13956995,   'charge': 1},
        'pi0':     {'geant_pid': 7,  'pdg_type': 111,        'mass': 0.1349768,    'charge': 0},
        'km':      {'geant_pid': 12, 'pdg_type': -321,       'mass': 0.493677,     'charge': -1},
        'kp':      {'geant_pid': 11, 'pdg_type': 321,        'mass': 0.493677,     'charge': 1},
        'eta':     {'geant_pid': 17, 'pdg_type': 221,        'mass': 0.547862,     'charge': 0},
        'proton':  {'geant_pid': 14, 'pdg_type': 2212,       'mass': 0.93827208,   'charge': 1},
        'pb208':   {'geant_pid': 111,'pdg_type': 1000822080, 'mass': 193.7507,     'charge': 82}
    }
    
    if particle_name not in particle_db:
        print(f"Error: Unknown particle '{particle_name}'")
        print(f"Known particles: {', '.join([k for k in particle_db.keys() if k not in ['photon', 'proton', 'pb208']])}")
        sys.exit(1)
    
    info = particle_db[particle_name].copy()
    info['name'] = particle_name
    return info


def ascii2hddm(final_state_string, target_type, input_file, output_file, 
               runNumber=0, vertex=None, use_doubles=True, include_recoil=False):
    """
    Convert ASCII generator output to HDDM format with arbitrary final state
    
    Args:
        final_state_string: String describing final state (e.g., "epem", "pi0pippim", "em")
        target_type: Target material - 'p' (proton), 'pb208' (lead-208)
        input_file: Input ASCII file path
        output_file: Output HDDM file path
        runNumber: Run number to embed in HDDM
        vertex: Vertex string 'vx vy zmin zmax'
        use_doubles: Add momentum_double for all particles (default: True)
        include_recoil: Write recoil nucleus to HDDM (default: False for pb208, True for proton)
    """
    
    # Parse final state
    particles = parse_final_state(final_state_string)
    n_particles = len(particles)
    
    # Get particle information
    particle_infos = [get_particle_info(p) for p in particles]
    
    # Validate target
    if target_type not in ['p', 'pb208']:
        print(f"Error: Invalid target_type '{target_type}'. Must be 'p' or 'pb208'.")
        sys.exit(1)
    
    # Target configuration
    target_name = 'proton' if target_type == 'p' else 'pb208'
    target_info = get_particle_info(target_name)
    
    # Determine if recoil should be written to HDDM
    # Default: include recoil for proton, exclude for pb208 (heavy nuclei never escape target)
    # User can override with include_recoil flag
    write_recoil = include_recoil if include_recoil else (target_type == 'p')
    
    # Calculate expected column count
    # Format: E_beam, [E1, p1x, p1y, p1z], [E2, p2x, p2y, p2z], ..., [E_recoil, prx, pry, prz]
    expected_cols = 1 + 4 * n_particles + 4  # beam + N*(E,px,py,pz) + recoil(E,px,py,pz)
    
    print(f"Final state: {final_state_string}")
    print(f"  Parsed as: {' + '.join(particles)} + recoil")
    if write_recoil:
        print(f"  HDDM products: {n_particles} particles + 1 recoil = {n_particles + 1} total")
    else:
        print(f"  HDDM products: {n_particles} particles (recoil not written to HDDM)")
    print(f"  Expected columns: {expected_cols}")
    print(f"    Format: E_beam, {', '.join([f'{p}(E,px,py,pz)' for p in particles])}, recoil(E,px,py,pz)")

    # Vertex configuration
    # Default: '0:0:0:0' tells downstream GlueX simulation (hdgeant4) to pull
    # vertex information from the beam photon and target specifications for the run.
    # Non-zero values override with fixed vertex position (vx, vy) and z-range (zmin, zmax).
    # The z-position will be randomly sampled uniformly between zmin and zmax for each event.
    vx, vy, vz_min, vz_max = 0.0, 0.0, 0.0, 0.0
    if vertex is not None:
        try:
            # Auto-detect delimiter: prefer colon (shell-safe), but support spaces for backward compatibility
            if ':' in vertex:
                vx, vy, vz_min, vz_max = map(float, vertex.split(':'))
            else:
                vx, vy, vz_min, vz_max = map(float, vertex.split())
        except Exception:
            print("Error parsing vertex. Format should be: 'vx:vy:zmin:zmax' (colon-separated, preferred) or 'vx vy zmin zmax' (space-separated)")
            sys.exit(1)

    # Open output stream
    fout = hddm_s.ostream(output_file)
    eventNo = 0
    skipped_lines = 0
    
    # Process input file
    for line_num, line in enumerate(open(input_file), 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        try:
            fixed_line = fix_exponent_format(line)
            fields = [float(f) for f in fixed_line.split(',')]
        except ValueError as e:
            print(f"Warning: Could not parse line {line_num}: {e}")
            skipped_lines += 1
            continue
        
        if len(fields) != expected_cols:
            if skipped_lines < 5:  # Only print first few errors
                print(f"Warning: Line {line_num} has {len(fields)} columns (expected {expected_cols}), skipping")
            skipped_lines += 1
            continue

        # Parse fields
        E_beam = fields[0]
        
        # Extract particle 4-momenta
        particle_4momenta = []
        for i in range(n_particles):
            start_idx = 1 + i * 4
            E, px, py, pz = fields[start_idx:start_idx + 4]
            particle_4momenta.append((E, px, py, pz))
        
        # Extract recoil 4-momentum
        recoil_start = 1 + n_particles * 4
        E_recoil, px_recoil, py_recoil, pz_recoil = fields[recoil_start:recoil_start + 4]

        # Create HDDM event structure
        rec = hddm_s.HDDM()
        eventNo += 1
        pev = rec.addPhysicsEvents()
        pev[0].runNo = runNumber
        pev[0].eventNo = eventNo
        
        rea = pev[0].addReactions()
        rea[0].type = 0
        rea[0].weight = 1

        # Beam photon
        bea = rea[0].addBeams()
        bea[0].type = get_particle_info('photon')['geant_pid']
        bmom = bea[0].addMomenta()
        bmom[0].E = float(np.float32(E_beam))
        bmom[0].px = 0
        bmom[0].py = 0
        bmom[0].pz = float(np.float32(E_beam))
        
        # Add double precision momentum for beam (default behavior)
        if use_doubles:
            bmom_double = bmom[0].addMomentum_doubles()
            bmom_double[0].E = E_beam
            bmom_double[0].px = 0.0
            bmom_double[0].py = 0.0
            bmom_double[0].pz = E_beam
        
        bpro = bea[0].addPropertiesList()
        bpro[0].charge = 0
        bpro[0].mass = 0

        # Target
        tar = rea[0].addTargets()
        tar[0].type = target_info['geant_pid']
        tmom = tar[0].addMomenta()
        tmom[0].E = target_info['mass']
        tmom[0].px = 0
        tmom[0].py = 0
        tmom[0].pz = 0
        tpol = tar[0].addPolarizations()
        tpol[0].Px = 0
        tpol[0].Py = 0
        tpol[0].Pz = 0
        tpro = tar[0].addPropertiesList()
        tpro[0].charge = target_info['charge']
        tpro[0].mass = target_info['mass']

        # Vertex
        vtx = rea[0].addVertices()
        origin = vtx[0].addOrigins()
        origin[0].vx = vx
        origin[0].vy = vy
        # If zmin == zmax, all events get the same z-value (useful for default '0 0 0 0')
        origin[0].vz = random.uniform(vz_min, vz_max)
        origin[0].t = 0.0
        
        # Create products: N particles + recoil (if writing recoil)
        n_products = n_particles + 1 if write_recoil else n_particles
        prod = vtx[0].addProducts(n_products)

        # Add all particles in order
        for i, (info, (E_particle, px_particle, py_particle, pz_particle)) in enumerate(zip(particle_infos, particle_4momenta)):
            prod[i].id = i + 1
            prod[i].pdgtype = info['pdg_type']
            prod[i].type = info['geant_pid']
            prod[i].mech = 909
            prod[i].parentid = 0
            prod[i].decayVertex = 0
            
            mom = prod[i].addMomenta()
            # PASS THROUGH values directly from generator
            # The generator has already enforced energy conservation in float32
            # DO NOT recalculate E - that would break energy conservation!
            mom[0].E = float(np.float32(E_particle))
            mom[0].px = float(np.float32(px_particle))
            mom[0].py = float(np.float32(py_particle))
            mom[0].pz = float(np.float32(pz_particle))
            
            # Add double precision momentum for ALL particles (default behavior)
            if use_doubles:
                mom_double = mom[0].addMomentum_doubles()
                mom_double[0].E = E_particle
                mom_double[0].px = px_particle
                mom_double[0].py = py_particle
                mom_double[0].pz = pz_particle
            
            # Add properties with correct mass
            props = prod[i].addPropertiesList()
            props[0].charge = info['charge']
            props[0].mass = info['mass']

        # Add recoil as final product (if requested)
        if write_recoil:
            recoil_idx = n_particles
            prod[recoil_idx].id = recoil_idx + 1
            prod[recoil_idx].pdgtype = target_info['pdg_type']
            prod[recoil_idx].type = target_info['geant_pid']
            prod[recoil_idx].mech = 909
            prod[recoil_idx].parentid = 0
            prod[recoil_idx].decayVertex = 0
            
            mom_recoil = prod[recoil_idx].addMomenta()
            # PASS THROUGH values directly from generator
            # The generator has already enforced energy conservation in float32
            # DO NOT recalculate E - that would break energy conservation!
            mom_recoil[0].E = float(np.float32(E_recoil))
            mom_recoil[0].px = float(np.float32(px_recoil))
            mom_recoil[0].py = float(np.float32(py_recoil))
            mom_recoil[0].pz = float(np.float32(pz_recoil))
            
            # Add double precision momentum for recoil (default behavior)
            if use_doubles:
                mom_recoil_double = mom_recoil[0].addMomentum_doubles()
                mom_recoil_double[0].E = E_recoil
                mom_recoil_double[0].px = px_recoil
                mom_recoil_double[0].py = py_recoil
                mom_recoil_double[0].pz = pz_recoil
            
            # Add properties with correct mass for recoil
            props_recoil = prod[recoil_idx].addPropertiesList()
            props_recoil[0].charge = target_info['charge']
            props_recoil[0].mass = target_info['mass']
        
        fout.write(rec)

    del fout
    
    print(f"\nConversion complete:")
    print(f"  Events written: {eventNo}")
    if skipped_lines > 0:
        print(f"  Lines skipped: {skipped_lines}")
    print(f"  Output file: {output_file}")
    print(f"  Products per event: {n_products} ({n_particles} particles{' + recoil' if write_recoil else ', no recoil'})")
    if use_doubles:
        print(f"  Precision: Double precision (momentum_double added)")
    else:
        print(f"  Precision: Single precision only")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ASCII to HDDM converter - specify any final state",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Particle Codes:
  em   = e-         (electron)
  ep   = e+         (positron)
  mum  = mu-        (muon)
  mup  = mu+        (antimuon)
  pim  = pi-        (negative pion)
  pip  = pi+        (positive pion)
  pi0  = pi0        (neutral pion)
  km   = K-         (negative kaon)
  kp   = K+         (positive kaon)
  eta  = eta        (eta meson)

Examples:
  %(prog)s epem p input.txt output.hddm              # e+ e- on proton (includes recoil proton)
  %(prog)s em p single_e.txt output.hddm             # single electron on proton
  %(prog)s mupmum pb208 muons.txt output.hddm        # mu+ mu- on lead (NO recoil by default)
  %(prog)s epem pb208 input.txt out.hddm --include-recoil  # Force recoil Pb208 in HDDM
  %(prog)s pi0 p pi0_events.txt output.hddm          # single pi0
  %(prog)s pi0pippim p three_pion.txt output.hddm    # pi0 pi+ pi- (3 pions)
  %(prog)s pippimpi0 p three_pion.txt output.hddm    # pi+ pi- pi0 (different order)
  %(prog)s ememem p triple_e.txt output.hddm         # 3 electrons
  %(prog)s epempippim p mixed.txt output.hddm        # e+ e- pi+ pi-

Recoil handling:
  - Proton target (p): Recoil proton written to HDDM by default
  - Lead target (pb208): Recoil Pb208 NOT written by default (saves ~10-15% file size)
  - Use --include-recoil to force writing recoil for any target

The script automatically determines:
  - Number of particles (from final state string)
  - Expected column count (1 + 4*N + 4)
  - Product ordering (as specified in final state)
        """)
    
    parser.add_argument("final_state", 
                        help="Final state particle string (e.g., 'epem', 'pi0', 'mupmumpippim')")
    parser.add_argument("target", 
                        choices=["p", "pb208"], 
                        help="Target: p (proton) or pb208 (lead-208)")
    parser.add_argument("input_file", 
                        help="Input ASCII file (comma-separated)")
    parser.add_argument("output_file", 
                        help="Output HDDM file")
    parser.add_argument("--run", 
                        type=int, default=0, 
                        help="Run number (default: 0)")
    parser.add_argument("--vertex", 
                        type=str, 
                        default="0:0:0:0",
                        help="Vertex position as 'vx:vy:zmin:zmax' (colon-separated, preferred) or 'vx vy zmin zmax' (space-separated). "
                             "Default '0:0:0:0' tells hdgeant4 to use run-specific beam/target geometry. "
                             "Non-zero values specify fixed vertex (vx,vy) and z-range for uniform sampling.")
    parser.add_argument("--no-doubles", 
                        action="store_true",
                        help="Do NOT add momentum_double (use single precision only). "
                             "Default is to add double precision for all momenta.")
    parser.add_argument("--include-recoil", 
                        action="store_true",
                        help="Write recoil nucleus to HDDM file. "
                             "Default: True for proton, False for pb208 (heavy recoils never escape target).")
    
    args = parser.parse_args()

    ascii2hddm(
        final_state_string=args.final_state,
        target_type=args.target,
        input_file=args.input_file,
        output_file=args.output_file,
        use_doubles=not args.no_doubles,
        include_recoil=args.include_recoil,
        runNumber=args.run,
        vertex=args.vertex
    )
