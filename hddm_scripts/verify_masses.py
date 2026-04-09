#!/usr/bin/env python3
"""
Verify masses in HDDM file by calculating from stored 4-momentum.
Uses momentum_double (double precision) if available, falls back to float32.

Usage: python3 verify_masses.py <hddm_file> [--events N]
"""
import hddm_s
import sys
import numpy as np

# Parse args
hddm_file = None
max_events = None  # None = all events
args = sys.argv[1:]
while args:
    a = args.pop(0)
    if a == '--events' and args:
        max_events = int(args.pop(0))
    elif not hddm_file:
        hddm_file = a

if not hddm_file:
    print("Usage: python3 verify_masses.py <hddm_file> [--events N]")
    sys.exit(1)

print(f"Checking masses in {hddm_file}...")
if max_events:
    print(f"(limiting to first {max_events} events)")
print()

electron_mass = 0.000510999  # GeV
muon_mass     = 0.105658389  # GeV
proton_mass   = 0.93827208   # GeV

em_masses = []
ep_masses = []
mup_masses = []
mum_masses = []
recoil_masses = []

using_doubles = None

# Read events
for event_num, rec in enumerate(hddm_s.istream(hddm_file)):
    if max_events is not None and event_num >= max_events:
        break
    physicsEvent = rec.getPhysicsEvents()[0]
    reaction = physicsEvent.getReactions()[0]
    vertex = reaction.getVertices()[0]
    products = vertex.getProducts()
    
    for product in products:
        pdg = product.pdgtype
        momentum = product.getMomenta()[0]

        # Prefer momentum_double over float32
        mom_doubles = momentum.getMomentum_doubles()
        if len(mom_doubles) > 0:
            md = mom_doubles[0]
            E  = md.E
            px = md.px
            py = md.py
            pz = md.pz
            if using_doubles is None:
                using_doubles = True
        else:
            E  = momentum.E
            px = momentum.px
            py = momentum.py
            pz = momentum.pz
            if using_doubles is None:
                using_doubles = False
        
        # Calculate mass: m² = E² - p²
        p_squared = px**2 + py**2 + pz**2
        m_squared = E**2 - p_squared
        
        if m_squared < 0:
            mass = -np.sqrt(-m_squared)
        else:
            mass = np.sqrt(m_squared)
        
        mass_mev = mass * 1000  # Convert to MeV
        
        if pdg == 11:   em_masses.append(mass_mev)
        elif pdg == -11:  ep_masses.append(mass_mev)
        elif pdg == 13:   mum_masses.append(mass_mev)
        elif pdg == -13:  mup_masses.append(mass_mev)
        elif pdg == 2212: recoil_masses.append(mass_mev)

# Statistics
n_events = event_num + 1 if 'event_num' in dir() else 0
print(f"Analyzed {n_events} events")
prec_src = 'momentum_double (double)' if using_doubles else 'standard momentum (float32) -- doubles not found!'
print(f"Precision source: {prec_src}")
print()

def print_stats(label, masses, expected_mev):
    if not masses:
        return
    arr = np.array(masses)
    errors = arr - expected_mev
    print(f"{label}:")
    print(f"  Expected: {expected_mev:.6f} MeV")
    print(f"  Mean:     {np.mean(arr):.6f} MeV  (error: {np.mean(errors):+.6f} MeV)")
    print(f"  Std dev:  {np.std(arr):.6f} MeV")
    print(f"  Error range: {np.min(errors):+.6f} to {np.max(errors):+.6f} MeV")
    print()

print_stats("Electron (e-, PDG 11)",  em_masses,  electron_mass * 1000)
print_stats("Positron (e+, PDG -11)", ep_masses,  electron_mass * 1000)
print_stats("Muon- (PDG 13)",         mum_masses, muon_mass     * 1000)
print_stats("Muon+ (PDG -13)",        mup_masses, muon_mass     * 1000)
print_stats("Recoil proton",          recoil_masses, proton_mass * 1000)
