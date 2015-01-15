#rocket-calc.py
"""

TODO: Everything needs to be reorganized to fit with the Click api
Script for Basic Calulations

Usage:
<insert documentation here>

"""
import sys
import getopt
import math

class constants:
	"""

	R = gas constant <units needed>
	y = heat capacity ratio or adiabatic index or ratio of specific heats
	g_c = earth gravitational constant <units needed>
	P_atm = atmospheric pressure, psi

	"""
	R = 65
	y = None
	g_c = 32.2
	P_atm = 14.7

	def __init__(self, fuel, oxidizer):
		if (fuel == "hydrocarbon" and oxidizer == "gaseous oxygen"):
			self.y = 1.2
		#Insert more fuel/oxidizer combinations here

class variables:
	"""

	r = O/F ratio (or mixture ratio)
	P_c = chamber pressure, psi sea level
	T_c = chamber gas temperature, °F/R
	L* = the chamber volume required for complete combustion (characteristic chamber length)
	F = thrust
	I_sp = specific impulse performance of fuel oxidizer combination as a result of chamber pressure and mixture ratio, seconds
	#TODO: C_d = orifice discharge coefficient (constant based on injector)
	#TODO: delta_P_injector = pressure drop across orifice, lb/ft2
	#TODO: rho_f = density of propellant, lb/ft3

	"""
	r, P_c, T_c, L, F, I_sp = None, None, None, None, None, None
	w_t, w_o, w_f = None, None, None

	def __init__(self, mixture_ratio, chamber_pressure, chamber_gas_temperature, characteristic_chamber_length, thrust, Isp):
		self.r = mixture_ratio
		self.P_c = chamber_pressure
		self.T_c = chamber_gas_temperature
		self.L = characteristic_chamber_length
		self.F = thrust
		self.I_sp = Isp
		#Derived From Given
		self.w_t = self.total_flow_rate()
		self.w_o = self.oxidizer_flow_rate()
		self.w_f = self.fuel_flow_rate()

	"""

	Design Equation:
		w_t = F/I_sp

	Variables:
		w_t = total flow rate, lb/sec
		F = thrust
		I_sp = specific impulse, seconds

	"""
	def total_flow_rate(self):
		return self.F/self.I_sp
	"""

	Design Equation:
		w_o = w_t*r / (r + 1)

	Variables:
		w_o = oxidizer flow rate, lb/sec
		w_t = total flow rate, lb/sec
		r = O/F ratio (or mixture ratio)
	"""
	def oxidizer_flow_rate(self):
		return (self.total_flow_rate()*self.r)/(self.r + 1)
	"""

	Design Equation:
		w_f = w_t / (r + 1)

	Variables:
		w_f = fuel flow rate, lb/sec
		w_t = total flow rate, lb/sec
		r = O/F ratio (or mixture ratio)
	"""
	def fuel_flow_rate(self):
		return (self.total_flow_rate())/(self.r + 1)

# class variables:
# 	"""

# 	Q = total heat transferred, Btu/sec
# 	q = average heat transfer rate of chamber, Btu/in2-sec
# 	A = heat transfer area, in2
# 	w_w = coolant flow rate, lb/sec
# 	c_p = specific heat of coolant, Btu/lb°F
# 	T = temperature of coolant leaving jacket, °F
# 	T_i = temperature of coolant entering jacket, °F

# 	"""
# 	Q, q, A, w_w, c_p, T, T_i = None, None, None, None, None, None, None

class nozzle:
	"""

	T_t = temperature of the gases at the nozzle throat <units needed>
	P_t = gas pressure at the nozzle throat
	A_t = nozzle throat cross-sectional area <units needed>
	D_t = nozzle throat diameter, in
	A_e = the nozzle exit cross-sectional area corresponding to the exit Mach number resulting from the choice of chamber pressure <units needed>
	D_e = nozzle exit diameter, in

	"""
	T_t, P_t, A_t, D_t, A_e, D_e = None, None, None, None, None, None

	def __init__(self, variables, constants):
		self.T_t = self.tempurature_of_gases_at_nozzle_throat(variables.T_c, constants)
		self.P_t = self.gas_pressure_at_nozzle_throat(variables.P_c, constants)
		self.A_t = self.nozzle_throat_cross_sectional_area(variables.w_t, self.P_t, self.T_t, constants)
		self.D_t = self.nozzle_throat_diameter(self.A_t)
		self.M_e = self.mach_number(variables.P_c, constants) #broken
		self.A_e = self.nozzle_exit_cross_sectional_area(self.A_t, self.M_e, constants) #broken
		self.D_e = self.nozzle_exit_diameter(self.A_e)

	"""
	
	Design Equation:
		A_t = (w_t/P_t) * sqrt((R*T_t)/(y*g_c))

	Variables:
		A_t = nozzle throat cross-sectional area <units needed>
		w_t = total weight flow (fuel + oxidizer) lb/sec
			w_o = (w_t * r)/(r + 1) lb/sec
			w_f = (w_t)/(r + 1) lb/sec
			w_t = w_o + w_f lb/sec
		P_t = gas pressure at the nozzle throat <units needed>
			<insert equations>
		R = gas constant <units needed>
		T_t = temperature of the gases at the nozzle throat <units needed>
			<insert equations>
		y = heat capacity ratio or adiabatic index or ratio of specific heats
		g_c = earth gravitational constant <units needed>

	"""
	def nozzle_throat_cross_sectional_area(self, w_t, P_t, T_t, constants):
		return (w_t/P_t)*math.sqrt((constants.R*T_t)/(constants.y*constants.g_c))
	"""
	
	Design Equation:
		T_t = T_c [1 / (1 + ((y - 1)/2))]

	Variables:
		T_t = temperature of the gases at the nozzle throat <units needed>
		T_c = the combustion chamber flame temperature in degrees Rankine (°R)
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def tempurature_of_gases_at_nozzle_throat(self, T_c, constants):
		return T_c*(1/(1 + ((constants.y - 1)/2)))
	"""

	Design Equation: 
		P_t = P_c * [ 1 + (y-1)/(2) ]^(-(y)/(y-1))

	Variables:
		P_t = gas pressure at the nozzle throat
		P_c = pressure in the combustion chamber <units needed>
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def gas_pressure_at_nozzle_throat(self, P_c, constants):
		return P_c*((1 + ((constants.y-1)/2))**(-(constants.y)/(constants.y-1)))
	"""

	Design Equation:
		M_e^2 = ((2)/(y-1))*[(P_c/P_atm)^((y-1)/(y)) - 1] 

	Variables:
		M_e = the ratio of the gas velocity to the local speed of sound (mach number) at the nozzle exit
		P_c = the pressure in the combustion chamber, psi
		P_atm = atmospheric pressure, psi
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def mach_number(self, P_c, constants):
		return math.sqrt(((2)/(constants.y-1))*(((P_c/constants.P_atm)**((constants.y-1)/(constants.y))) - 1)) 
	"""

	Design Equation:
		A_e = (A_t/M_e)*[(1 + ((y - 1)/(2))*M_e^2)/((y + 1)/(2))]^((y+1)/(2*(y-1)))

	Variables:
		A_e = the nozzle exit cross-sectional area corresponding to the exit Mach number resulting from the choice of chamber pressure, in^2
		A_t = nozzle throat cross-sectional area, in^2
		M_e = the ratio of the gas velocity to the local speed of sound (mach number) at the nozzle exit
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def nozzle_exit_cross_sectional_area(self, A_t, M_e, constants):
		return (A_t/M_e)*((1 + ((constants.y - 1)/(2))*M_e**2)/((constants.y + 1)/(2)))**((constants.y+1)/(2*(constants.y-1)))
	"""

	Design Equation: 
		D_t = sqrt((4*A_t)/pi)

	Variables:
		D_t = nozzle throat diameter, in
		A_t = nozzle throat cross-sectional area, in^2

	"""
	def nozzle_throat_diameter(self, A_t):
		return math.sqrt((4*A_t)/(3.14))
	"""

	Design Equation:
		D_e = sqrt((4*A_e)/pi)

	Variables:
		D_e = nozzle exit diameter, in
		A_e = the nozzle exit cross-sectional area corresponding to the exit Mach number resulting from the choice of chamber pressure, in^2

	"""
	def nozzle_exit_diameter(self, A_e):
		return math.sqrt((4*A_e)/(3.14))

class combustion_chamber:
	"""

	V_c = the volume of the combustion chamber, in^3
	D_c = chamber diameter (we assume this is five times the nozzle throat diameter), in
	A_c = the combustion chamber cross-sectional area, in^2
	L_c = chamber length, in
	t_w = chamber wall thickness, in
	S = allowable working stress on the combustion chamber wall (psi) based on material, psi
	material = material combustion chamber is made of, only copper supported currently

	"""
	material = None
	V_c, D_c, A_c, L_c, S, t_w = None, None, None, None, None, None

	def __init__(self, variables, constants, nozzle, chamber_material):
		self.V_c = self.combustion_chamber_volume(nozzle.A_t, variables.L)
		self.D_c = nozzle.D_t*5
		self.A_c = self.combustion_chamber_cross_sectional_area(self.D_c)
		self.L_c = self.chamber_length(self.V_c, self.D_c, self.A_c)
		if (chamber_material == "copper"):
			self.S = 8000
		self.t_w = self.combustion_chamber_thickness(variables.P_c, self.D_c, self.S)
		self.material = chamber_material 

	"""

	Design Equation:
		L* = V_c/A_t

	Variables:
		L* = the chamber volume required for complete combustion (characteristic chamber length), in
		V_c = the chamber volume (including the converging section of the nozzle), in^3
		A_t = the nozzle throat area, in^2

	"""
	def characteristic_chamber_length(self, V_c, A_t):
		return (V_c/A_t)
	"""

	Design Equation:
		L_c = V_c/(D_c * A_c)

	Variables:
		L_c = chamber length, in
		D_c = chamber diameter (we assume this is five times the nozzle throat diameter), in
		A_c = the combustion chamber cross-sectional area, in^2

	"""
	def chamber_length(self, V_c, D_c, A_c):
		return ((V_c)/(D_c*A_c))
	"""

	Design Equation:
		A_c = (pi*D_c^2)/4

	Variables:
		A_c = the combustion chamber cross-sectional area, in^2
		D_c = the combustion chamber diameter, in

	"""
	def combustion_chamber_cross_sectional_area(self, D_c):
		return (3.14*D_c**2)/4
	"""
	
	Design Equation:
		V_c = A_c * L_c + convergent volume

	Variables:
		V_c = the volume of the combustion chamber, in^3
		A_c = the combustion chamber cross-sectional area, in^2
		L_c = length of the combustion chamber, in

	Note: If convergent volume is left out it approximates it to 1/10th the volume of the cylindrical portion chamber

	"""
	def combustion_chamber_volume(self, A_c, L_c, convergent_volume):
		return A_c*L_c + convergent_volume
	"""

	Design Equation:
		V_c = L* A_t

	Variables:
		V_c = the volume of the combustion chamber, in^3
		A_t = the nozzle throat area, in^2
		L* = the chamber volume required for complete combustion (characteristic chamber length), in

	"""
	def combustion_chamber_volume(self, A_t, L):
		return L*A_t
	"""

	Design Equation:
		S = P*D/2*t_w

	Variables:
		S = working stress on the combustion chamber wall
		P = the pressure in the combustion chamber (neglecting the effect of coolant pressure on the outside of the shell), psi
		D = the mean diameter of the cylinder, in
		t_w = the thickness of the cylinder wall, in

	Note: Chamber pressure is usually a controlled aspect of the design 

	"""
	def combustion_chamber_stress(self, P, D, t_w):
		return ((P*D)/(t_w*2))
	"""

	Design Equation:
		t_w = (P*D)/(2*S)

	Variables:
		t_w = the thickness of the cylinder wall
		P = the pressure in the combustion chamber (neglecting the effect of coolant pressure on the outside of the shell), psi
		D = the mean diameter of the cylinder, in
		S = the allowable working stress (psi) based on the material of the chamber, psi

	Note: This is the minimum thickness since welding factors and design considerations (such as O-ring grooves, etc.) will 
	usually require walls thicker than those indication by the stress equation

	"""
	def combustion_chamber_thickness(self, P, D, S):
		return ((P*D)/(S*2))

class engine_cooling: #This may require some organizational rethinking. Maybe coolant should be selected at the beginning. 
	"""

	A = heat transfer area, in^2
	q = average heat transfer rate for coolant, Btu/in^2-sec
	Q = total heat transferred, Btu/sec
	w_cool = coolant flow rate, lb/sec
	v_cool = flow velocity of the coolant, ft/sec
	D_1 = the outer diameter of the combustion chamber, in
	D_2 = the inner diameter of the outer jacket, in
	delta_D = annular flow passage width, in
	rho_cool = density of coolant, lb/ft^3
	delta_T = allowable/desired rise in temperature of coolant (defaults to 40), °F

	"""
	A, q, Q, w_cool, v_cool, D_1, D_2, delta_D, rho_cool, delta_T = None, None, None, None, None, None, None, None, None, None

	def __init__(self, variables, constants, nozzle, combustion_chamber, coolant, coolant_tempurature_rise=40, v_cool=30):
		if(coolant == "water"):
			self.rho_cool = 62.4
		self.coolant_tempurature_rise = coolant_tempurature_rise

		self.A = self.heat_transfer_area(combustion_chamber.D_c, combustion_chamber.t_w, combustion_chamber.L_c, 0)*1.1
		if(combustion_chamber.material == "copper"):
			self.q = 3
		self.Q = self.total_heat_transfer_from_chamber_to_coolant(self.q, self.A)
		self.w_cool = self.coolant_flow_rate(self.Q, self.coolant_tempurature_rise)
		self.v_cool = v_cool
		self.D_1 = self.combustion_chamber_outer_diameter(combustion_chamber.D_c, combustion_chamber.t_w) #1.3875
		self.D_2 = self.cooling_jacket_inner_diameter(self.w_cool, self.v_cool, self.rho_cool, self.D_1) #1.475
		self.delta_D = self.D_2 - self.D_1 

	"""

	Design Equation:
		Q = q*A = w_w * c_p (T - T_i)

	Variables:
		Q = total heat transferred, Btu/sec
		q = average heat transfer rate of chamber, Btu/in^2-sec
		A = heat transfer area, in^2
		w_w = coolant flow rate, lb/sec
		c_p = specific heat of coolant, Btu/lb°F
		T = temperature of coolant leaving jacket, °F
		T_i = temperature of coolant entering jacket, °F

	"""
	def total_heat_transfer_from_chamber_to_coolant(self, w_w, c_p, T, T_i):
		return ((w_w*c_p)*(T-T_i))
	"""

	Design Equation:
		Q = q*A

	Variables:
		Q = total heat transferred, Btu/sec
		q = average heat transfer rate of chamber, Btu/in^2-sec
		A = heat transfer area, in^2

	"""
	def total_heat_transfer_from_chamber_to_coolant(self, q, A):
		return q*A
	"""

	Design Equation:
		A = pi (D_c + 2*t_w)(L_c) + area of nozzle cone

	Variables:
		A = heat transfer area, in^2
		D_c = combustion chamber diameter, in
		t_w = combustion chamber wall thickness, in
		L_c = combustion chamber length, in

	Note: the area of the nozzle cone up to the throat can be assumed to be about 10 percent of the chamber surface area
	if it is not known 

	"""
	def heat_transfer_area(self, D_c, t_w, L_c, area_of_nozzle_cone):
		return 3.14*(D_c + 2*t_w)*L_c + area_of_nozzle_cone
	"""

	Design Equation:
		w_cool = Q/delta_T

	Variables:
		w_cool = coolant flow rate, lb/sec
		Q = total heat transferred, Btu/sec
		delta_T = desired temperature rise of the coolant, °F

	"""
	def coolant_flow_rate(self, Q, delta_T):
		return (Q/delta_T)
	"""

	Design Equation:
		v_cool = w_cool/rho_cool*A

	Variables:
		v_cool = flow velocity of the coolant, ft/sec
		w_cool = coolant flow rate, lb/sec
		A = heat transfer area, in^2

	"""
	def coolant_velocity(self, w_cool, rho_cool, A):
		return w_cool/(rho_cool*A)
	"""

	Design Equation:
		D_1 = D_c + 2*t_w

	Variables:
		D_1 = the outer diameter of the combustion chamber, in
		D_c = combustion chamber diameter, in
		t_w = combustion chamber wall thickness, in

	"""
	def combustion_chamber_outer_diameter(self, D_c, t_w):
		return (D_c + 2*t_w)
	"""

	Design Equation:
		D_2 = sqrt((4*w_cool/v_cool*rho_cool*pi) + (D_1^2))

	Variables:
		D_2 = the inner diameter of the outer jacket, in
		w_cool = coolant flow rate, lb/sec
		v_cool = flow velocity of the coolant, ft/sec
		rho_cool = density of coolant, lb/ft^3
		D_1 = the outer diameter of the combustion chamber, in

	"""
	def cooling_jacket_inner_diameter(self, w_cool, v_cool, rho_cool, D_1):
		return math.sqrt(((4*w_cool)/(v_cool*rho_cool*3.14)) + (D_1/12)**2)*12 #Convert from ft to inches

class injector: #Not Complete
	"""

	A_impinging = flow area for fuel injection if an impinging jet injector is chosen, ft^2
	A_spray = flow area for fuel injection if a spray nozzle injector is chosen, ft^2
	n_f = number of fuel injection holes
	D_f = diameter of fuel injection holes based on n
	n_o = number of oxidizer injection holes
	D_o = diameter of oxidizer injection holes based on n
	C_d = orifice discharge coefficient (constant based on injector)
	v_o = oxidizer injection velocity, ft/sec
	v_f = fuel injection velocity, ft/sec

	"""

	A_impinging, A_spray, n_f, D_f, n_o, D_o, C_d, v_o, v_f = None, None, None, None, None, None, None, None, None

	"""

	Design Equation:
		w = C_d * A sqrt(2 * g_c * rho_f * delta_P)

	Variables:
		w = propellant flow rate, lb/sec
		C_d = orifice discharge coefficient
		A = area of orifice, ft^2
		g_c = gravitational constant, 32.2 ft/sec^2
		rho_f = density of propellant, lb/ft^3
		delta_P = pressure drop across orifice, lb/ft^2

	"""
	def propellant_flow_rate(self, A, rho_f, delta_P, C_d=0.7):
		return (C_d*A)*math.sqrt(2*32.2*rho*delta_P)
	"""
	
	Design Equation:
		A_impinging = (w_f)/((C_d) * (2*g_c * rho_f * delta_P_injector)^(1/2))

	Variables:
		A_impinging = flow area for fuel injection if an impinging jet injector is chosen, ft^2
		w_f = propellant flow rate, lb/sec
		C_d = orifice discharge coefficient
		g_c = gravitational constant, 32.2 ft/sec2
		rho_f = density of propellant, lb/ft3
		delta_P_injector = pressure drop across orifice, lb/ft2

	"""
	def injector_flow_area_impinging(self, w_f, rho_f, delta_P_injector, C_d=0.7):
		return (w_f)/((C_d)*math.sqrt(2*g_c*rho_f*delta_P_injector))
	"""

	Design Equation:
		A_spray = w_o/rho_o*v_o

	Variables:
		A_spray = flow area for fuel injection if a spray nozzle injector is chosen, ft^2
		w_o = oxidizer flow rate, lb/sec
		rho_o = density of oxidizer inside injector, lb/ft3
		v_o = oxidizer injection velocity, ft/sec

	Related Equation: The Ideal Gas Law at constant tempurature is used to find rho_f
		rho_o_inside = rho_o_outside * ( (delta_P_injector + P_c)/(P_atm) )

	Variables:
		delta_P_injector = pressure drop across orifice, lb/ft2
		P_c = chamber pressure, psi
		P_atm = atmospheric pressure, psi

	"""
	def injector_flow_area_spray(self, w_o, delta_P_injector, P_c, v_o, rho_o=0.083, P_atm=14.7):
		return (w_o)/((rho_o*((delta_P_injector+P_c)/(P_atm)))*(v_o))

	"""

	Design Equation:
		v = C_d * sqrt(2*g_c * (delta_P/rho))

	Variables:
		v = the injection velocity, or velocity of the liquid stream issuing from the orifice, ft/sec
		C_d = orifice discharge coefficient
		g_c = gravitational constant, 32.2 ft/sec^2
		delta_P = pressure drop across orifice, lb/ft^2
		rho_f = density of propellant, lb/ft^3

	"""
	def injection_velocity(self, delta_P, rho, C_d=0.6):
		return C_d*math.sqrt((2*g_c)*(delta_P/rho))

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def main():
    print(bcolors.HEADER + "===========Begin Report===========" + bcolors.ENDC)
    """
    
    Known Values:
    	mixture_ratio = 2.5
    	chamber_pressure = 300 psi
    	chamber_gas_temperature = 6202 R
    	characteristic_chamber_length = 60 in
    	thrust = 20lb
    	Isp = 260 sec

    """
    var = variables(2.5, 300, 6202, 60, 20, 260)  
    con = constants("hydrocarbon", "gaseous oxygen")

    print(bcolors.OKGREEN + "===========Flow Rates===========" + bcolors.ENDC)
    print("total flow rate: " + str(var.w_t) + " lb/sec")
    print("fuel flow rate: " + str(var.w_f) + " lb/sec")
    print("oxidizer flow rate: " + str(var.w_o) + " lb/sec")

    print(bcolors.WARNING + "===========Nozzle===========" + bcolors.ENDC)
    noz = nozzle(var, con)
    print("tempurature of gases at nozzle throat: " + str(noz.T_t) + " °R")
    print("gas presure at nozzle throat: " + str(noz.P_t) + " psi")
    print("nozzle throat cross-sectional area: " + str(noz.A_t) + " in^2")
    print("nozzle throat diameter: " + str(noz.D_t) + " in")
    print("mach number: " + str(noz.M_e) + " (unitless?)")
    print("nozzle exit cross-sectional area: " + str(noz.A_e) + " in^2")
    print("nozzle exit diameter: " + str(noz.D_e) + " in")

    print(bcolors.FAIL + "===========Combustion Chamber===========" + bcolors.ENDC)
    com = combustion_chamber(var, con, noz, "copper")
    print("combustion chamber volume: " + str(com.V_c) + " in^3")
    print("combustion chamber diameter: " + str(com.D_c) + " in")
    print("combustion chamber area: " + str(com.A_c) + " in^2")
    print("combustion chamber length: " + str(com.L_c) + " in")
    print("combustion chamber wall thickness: " + str(com.t_w) + " in")

    print(bcolors.OKBLUE + "===========Engine Cooling===========" + bcolors.ENDC) #This looks like it has some problems
    eng_cool = engine_cooling(var, con, noz, com, "water")
    print("heat transfer area: " + str(eng_cool.A))
    print("average heat transfer rate for coolant: " + str(eng_cool.q))
    print("total heat transferred: " + str(eng_cool.Q) + " Btu/sec")
    print("coolant flow rate: " + str(eng_cool.w_cool))
    print("flow velocity of coolant: " + str(eng_cool.v_cool))
    print("outer diameter of the combustion chamber: " + str(eng_cool.D_1))
    print("inner diameter of cooling jacket: " + str(eng_cool.D_2))
    print("annular flow passage width: " + str(eng_cool.delta_D))

    print(bcolors.WARNING + "===========Injector===========" + bcolors.ENDC)


if __name__ == "__main__":
    main()