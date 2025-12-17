#!/usr/bin/env python3
"""
Verify masses in HDDM file by calculating from stored 4-momentum
"""
import hddm_s
import sys
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python3 verify_masses.py <hddm_file>")
    sys.exit(1)

hddm_file = sys.argv[1]

print(f"Checking masses in {hddm_file}...")
print()

electron_mass = 0.000510999  # GeV
proton_mass = 0.93827208     # GeV

em_masses = []
ep_masses = []
recoil_masses = []

# Read events
for rec in hddm_s.istream(hddm_file):
    physicsEvent = rec.getPhysicsEvents()[0]
    reaction = physicsEvent.getReactions()[0]
    vertex = reaction.getVertices()[0]
    products = vertex.getProducts()
    
    for product in products:
        pdg = product.pdgtype
        momentum = product.getMomenta()[0]
        
        E = momentum.E
        px = momentum.px
        py = momentum.py
        pz = momentum.pz
        
        # Calculate mass: m² = E² - p²
        p_squared = px**2 + py**2 + pz**2
        m_squared = E**2 - p_squared
        
        if m_squared < 0:
            mass = -np.sqrt(-m_squared)
        else:
            mass = np.sqrt(m_squared)
        
        mass_mev = mass * 1000  # Convert to MeV
        
        if pdg == 11:  # electron
            em_masses.append(mass_mev)
        elif pdg == -11:  # positron
            ep_masses.append(mass_mev)
        elif pdg == 2212:  # proton
            recoil_masses.append(mass_mev)

# Statistics
print(f"Analyzed {len(em_masses)} events")
print()

if em_masses:
    em_arr = np.array(em_masses)
    print(f"Electron (e-) masses:")
    print(f"  Expected: {electron_mass * 1000:.6f} MeV")
    print(f"  Mean:     {np.mean(em_arr):.6f} MeV")
    print(f"  Std dev:  {np.std(em_arr):.6f} MeV")
    print(f"  Min:      {np.min(em_arr):.6f} MeV")
    print(f"  Max:      {np.max(em_arr):.6f} MeV")
    errors = em_arr - electron_mass * 1000
    print(f"  Error range: {np.min(errors):.6f} to {np.max(errors):.6f} MeV")
    print()

if ep_masses:
    ep_arr = np.array(ep_masses)
    print(f"Positron (e+) masses:")
    print(f"  Expected: {electron_mass * 1000:.6f} MeV")
    print(f"  Mean:     {np.mean(ep_arr):.6f} MeV")
    print(f"  Std dev:  {np.std(ep_arr):.6f} MeV")
    print(f"  Min:      {np.min(ep_arr):.6f} MeV")
    print(f"  Max:      {np.max(ep_arr):.6f} MeV")
    errors = ep_arr - electron_mass * 1000
    print(f"  Error range: {np.min(errors):.6f} to {np.max(errors):.6f} MeV")
    print()

if recoil_masses:
    rec_arr = np.array(recoil_masses)
    print(f"Recoil proton masses:")
    print(f"  Expected: {proton_mass * 1000:.6f} MeV")
    print(f"  Mean:     {np.mean(rec_arr):.6f} MeV")
    print(f"  Std dev:  {np.std(rec_arr):.6f} MeV")
    print(f"  Min:      {np.min(rec_arr):.6f} MeV")
    print(f"  Max:      {np.max(rec_arr):.6f} MeV")
    errors = rec_arr - proton_mass * 1000
    print(f"  Error range: {np.min(errors):.6f} to {np.max(errors):.6f} MeV")
