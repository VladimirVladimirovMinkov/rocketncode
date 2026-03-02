#import softposit as sp


R = 8.31446261815324 #exakt

def solve(P=None, V=None, n=None, T=None, R=8.314):
    """Solves for the missing variable in PV = nRT."""
    if P is None: return (n * R * T) / V  # P = nRT / V
    if V is None: return (n * R * T) / P  # V = nRT / P
    if n is None: return (P * V) / (R * T) # n = PV / RT
    if T is None: return (P * V) / (n * R) # T = PV / nR

class gas:
    name = "air, dry"
    Pg = 101325 # sea level pressure
    Vg = 0.022414 # 1 mole at 1 atmo, m^3/1000 = L
    ng = 1.0 # number of moles
    Tk = 273.15 # 0 c
    mm = 28.96 # molar mass aka g/molc
    gamma = [28.11, 0.1967, 4.802, -1.966]
    
    def __init__(self, name = "air, dry", mm = 28.967, ng=1.0, Vg=0.022414, Tk=273.15, Pg=101325, cg = [28.11, 0.1967, 4.802, -1.966]):
        self.name = name
        self.Pg = Pg
        self.Vg = Vg
        self.ng = ng
        self.Tk = Tk
        self.mm = mm
        self.cg = cg

    #returns mass
    def mass(self):
        return self.ng*self.mm

    def dVg(self, Vg, process="thm"):
        process = process.lower()
        if process == "ithm":
            # P must change
            self.Pg = (self.ng * R * self.Tk) / Vg
            
        elif process == "ibar":
            # T must change
            self.Tk = (self.Pg * Vg) / (self.ng * R)
            
        elif process == "abat":
            # Both P and T change (gamma for air ~1.4)
            t = self.Tk/1000
            c = self.cg[0]+(self.cg[1]*t)+(self.cg[2]*t*t)+(self.cg[3]*t*t*t)
            gamma = c/(c-R)
            vr = self.Vg / Vg
            self.Pg = self.Pg * (vr ** gamma)
            self.Tk = self.Tk * (vr ** (gamma - 1)) # Use Adiabatic T-V relationship

        elif process == "cut":
            self.ng /= self.Vg/Vg

        self.Vg = Vg

    def dPg(self, Pg, process="ithm"):
        process = process.lower()
        if self.Pg == Pg: return
        pr = Pg / self.Pg  # Pressure ratio
        
        if process == "ithm":
            # P changes, V must compensate to keep T constant
            self.Vg = (self.ng * R * self.Tk) / Pg
            
        elif process == "abat":
            t = self.Tk / 1000
            cp = self.cg[0]+(self.cg[1]*t)+(self.cg[2]*t**2)+(self.cg[3]*t**3)
            gamma = cp / (cp - R)
            # V and T both change
            self.Vg = self.Vg * ( (1/pr) ** (1/gamma) )
            self.Tk = self.Tk * ( pr ** ((gamma - 1) / gamma) )

        self.Pg = Pg

    def dng(self, ng, process="ithm"):
        process = process.lower()
        if self.ng == ng: return
        
        if process == "ithm":
            # Volume is fixed (tank), P must change
            self.Pg = (ng * R * self.Tk) / self.Vg
            
        elif process == "ibar":
            # Pressure is fixed (balloon), V must change
            self.Vg = (ng * R * self.Tk) / self.Pg

        self.ng = ng

    def dTk(self, Tk, process="ibar"):
        process = process.lower()
        if self.Tk == Tk: return
        tr = Tk / self.Tk  # Temp ratio
        
        if process == "ibar":
            # P is fixed, V must change (Charles's Law)
            self.Vg = (self.ng * R * Tk) / self.Pg
            
        elif process == "abat":
            t = self.Tk / 1000
            cp = self.cg[0]+(self.cg[1]*t)+(self.cg[2]*t**2)+(self.cg[3]*t**3)
            gamma = cp / (cp - R)
            # V and P both change
            self.Vg = self.Vg * ( (1/tr) ** (1 / (gamma - 1)) )
            self.Pg = self.Pg * ( tr ** (gamma / (gamma - 1)) )

        self.Tk = Tk


    def print(self):
        print(self.name+ ":")
        print(f"Pressure:    Pa:  {self.Pg:<8.8g} Atmo: {self.Pg/101325:<8.8g} Bar: {self.Pg/100000:<8.8g} Tor: {self.Pg/133.322:<8.8g} Psi: {self.Pg/6894.76:<8.8g}")
        print(f"Moles:       n:   {self.ng:<8.8g} Kmol: {self.ng/1000:<8.8g}")
        print(f"Mass:        g:   {self.mass():<8.8g}   Kg: {self.mass()/1000:<8.8g}")
        print(f"Volume:      m^3: {self.Vg:<8.8g}    L: {self.Vg*1000:<8.8g}")
        print(f"Temperature: K:   {self.Tk:<8.8g}   °C: {self.Tk-273.15:<8.8g}  °F: {(self.Tk-273.15)*1.8+32:<8.8g}")
        print()

#P1V1/N1T1 = P2V2/N2T2

# Initialize the gases with your specific coefficients
air = gas()
co2 = gas(name="CO2", mm=44.009, cg=[24.99, 55.18, -33.69, 7.948])
n2o = gas(name="N2O", mm=44.013, cg=[26.85, 50.81, -31.42, 7.76])

gases = [air, co2, n2o]

# 1. TEST VOLUME (dVg) - ADIABATIC
print("--- TESTING dVg (Volume) | Process: abat ---")
for g in gases:
    g.print()
    g.dVg(g.Vg/2, "abat") # Half
    g.print()
    g.dVg(g.Vg*2, "abat") # Double (Back to start)
    g.print()
    print("-" * 30)

# 2. TEST PRESSURE (dPg) - ADIABATIC
print("\n--- TESTING dPg (Pressure) | Process: abat ---")
for g in gases:
    g.print()
    g.dPg(g.Pg/2, "abat") # Half
    g.print()
    g.dPg(g.Pg*2, "abat") # Double
    g.print()
    print("-" * 30)

# 3. TEST MOLES (dng) - ISOTHERMAL (Fixed Tank)
print("\n--- TESTING dng (Moles) | Process: ithm ---")
for g in gases:
    g.print()
    g.dng(g.ng/2, "ithm") # Half
    g.print()
    g.dng(g.ng*2, "ithm") # Double
    g.print()
    print("-" * 30)

# 4. TEST TEMPERATURE (dTk) - ISOBARIC (Expanding Volume)
print("\n--- TESTING dTk (Temperature) | Process: ibar ---")
for g in gases:
    g.print()
    g.dTk(g.Tk/2, "ibar") # Half
    g.print()
    g.dTk(g.Tk*2, "ibar") # Double
    g.print()
    print("-" * 30)

# 5. TEST THE "CUT" (Special Volume logic)
print("\n--- TESTING THE 'CUT' (V/2, n/2) ---")
for g in gases:
    g.print()
    g.dVg(g.Vg/2, "cut")
    g.print()
    # Note: "Cut" isn't naturally reversible by doubling V 
    # unless you define a "Join" process, but we can reset:
    g.dVg(g.Vg*2, "cut") 
    g.print()
    print("-" * 30)


def test_stepwise(gas_obj, target_V, process="ibar", steps=100):
    # Store initial state to compare later
    initial_n = gas_obj.ng
    
    # Calculate the step size
    total_delta_V = target_V - gas_obj.Vg
    step_size = total_delta_V / steps
    
    print(f"--- Step-wise Test: {gas_obj.name} | {process.upper()} | {steps} steps ---")
    gas_obj.print()
    
    # Perform the steps
    for i in range(steps):
        new_v = gas_obj.Vg + step_size
        gas_obj.dVg(new_v, process)
        
    print(f"Final state after {steps} steps:")
    gas_obj.print()
    
    # THE CROSS-CHECK: Using your solve() function
    # PV = nRT -> P should equal nRT/V
    expected_P = solve(V=gas_obj.Vg, n=gas_obj.ng, T=gas_obj.Tk, R=8.3144626)
    
    error_P = abs(gas_obj.Pg - expected_P)
    
    print(f"Ideal Gas Law Check (solve):")
    print(f"  Calculated P: {gas_obj.Pg:10.2f} Pa")
    print(f"  Expected P:   {expected_P:10.2f} Pa")
    print(f"  Accuracy Error: {error_P:10.6f} Pa")
    print("-" * 50)

# --- Usage Example ---
# Resetting air for a clean test
air = gas("Air Step Test")

# Compressing to 50% volume in 1 step vs 100 steps
test_stepwise(air, air.Vg/2, process="abat", steps=1)

# Reset and try with 100 steps (re-initialise gas to be sure)
air_precise = gas("Air High Precision")
test_stepwise(air_precise, air_precise.Vg/2, process="abat", steps=100)

