#rocket-calc.py
"""

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
	P_atm = atmospheric pressure (psi)

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
	P_c = chamber pressure (psi sea level)
	T_c = chamber gas temperature (°F/R)
	L* = the chamber volume required for complete combustion (characteristic chamber length)
	F = thrust
	I_sp = specific impulse performance of fuel oxidizer combination as a result of chamber pressure and mixture ratio

	"""
	r, P_c, T_c, L, F, I_sp = None, None, None, None, None, None

	def __init__(self, mixture_ratio, chamber_pressure, chamber_gas_temperature, characteristic_chamber_length, thrust, Isp):
		self.r = mixture_ratio
		self.P_c = chamber_pressure
		self.T_c = chamber_gas_temperature
		self.L = characteristic_chamber_length
		self.F = thrust
		self.I_sp = Isp 

	"""

	Design Equation:
		w_t = F/I_sp

	Variables:
		w_t = total flow rate
		F = thrust
		I_sp = specific impulse

	"""
	def total_flow_rate(self):
		return self.F/self.I_sp
	"""

	Design Equation:
		w_o = w_t*r / (r + 1)

	Variables:
		w_o = oxidizer flow rate
		w_t = total flow rate
		r = O/F ratio (or mixture ratio)
	"""
	def oxidizer_flow_rate(self):
		return (self.total_flow_rate()*self.r)/(self.r + 1)
	"""

	Design Equation:
		w_f = w_t / (r + 1)

	Variables:
		w_f = fuel flow rate
		w_t = total flow rate
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
		P_c = the pressure in the combustion chamber
		P_atm = atmospheric pressure (psi)
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def mach_number(self, P_c, constants):
		return ((2)/(constants.y-1))*(((P_c/constants.P_atm)**((constants.y-1)/(constants.y))) - 1) 
	"""

	Design Equation:
		A_e = (A_t/M_e)*[(1 + ((y - 1)/(2))*M_e^2)/((y + 1)/(2))]^((y+1)/(2*(y-1)))

	Variables:
		A_e = the nozzle exit cross-sectional area corresponding to the exit Mach number resulting from the choice of chamber pressure <units needed>
		A_t = nozzle throat cross-sectional area <units needed>
		M_e = the ratio of the gas velocity to the local speed of sound (mach number) at the nozzle exit
		y = heat capacity ratio or adiabatic index or ratio of specific heats

	"""
	def nozzle_exit_cross_sectional_area(self, A_t, M_e, constants):
		return (A_t/M_e)*((1 + ((constants.y - 1)/(2))*M_e**2)/((constants.y + 1)/(2)))**((constants.y+1)/(2*(constants.y-1)))

class combustion_chamber:
	"""

	Design Equation:
		L* = V_c/A_t

	Variables:
		L* = the chamber volume required for complete combustion (characteristic chamber length)
		V_c = the chamber volume (including the converging section of the nozzle), in cubic inches
		A_t = the nozzle throat area (in^2)

	"""
	def characteristic_chamber_length(self, V_c, A_t):
		return (V_c/A_t)
	"""

	Design Equation:
		A_c = (pi*D_c^2)/4

	Variables:
		A_c = the combustion chamber cross-sectional area
		D_c = the combustion chamber diameter

	"""
	def combustion_chamber_cross_sectional_area(self, D_c):
		return (3.14*D_c**2)/4
	"""
	
	Design Equation:
		V_c = A_c * L_c + convergent volume

	Variables:
		V_c = the volume of the combustion chamber
		A_c = the combustion chamber cross-sectional area
		L_c = length of the combustion chamber

	Note: If convergent volume is left out it approximates it to 1/10th the volume of the cylindrical portion chamber

	"""
	def combustion_chamber_volume(self, A_c, L_c, convergent_volume):
		return A_c* L_c + convergent_volume
	"""

	Design Equation:
		S = P*D/2*t_w

	Variables:
		S = working stress on the combustion chamber wall
		P = the pressure in the combustion chamber (neglecting the effect of coolant pressure on the outside of the shell)
		D = the mean diameter of the cylinder
		t_w = the thickness of the cylinder wall

	Note: Chamber pressure is usually a controlled aspect of the design 

	"""
	def combustion_chamber_stress(self, P, D, t_w):
		return ((P*D)/(t_w*2))
	"""

	Design Equation:
		t_w = (P*D)/(2*S)

	Variables:
		t_w = the thickness of the cylinder wall
		P = the pressure in the combustion chamber (neglecting the effect of coolant pressure on the outside of the shell)
		D = the mean diameter of the cylinder
		S = the allowable working stress (psi) based on the material of the chamber

	Note: This is the minimum thickness since welding factors and design considerations (such as O-ring grooves, etc.) will 
	usually require walls thicker than those indication by the stress equation

	"""
	def combustion_chamber_thickness(self, P, D, S):
		return ((P*D)/(S*2))

class engine_cooling:
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

class injector:
	"""

	Design Equation:
		w = C_d * A sqrt(2 * g_c * rho * delta_P)

	Variables:
		w = propellant flow rate, lb/sec
		C_d = orifice discharge coefficient
		A = area of orifice, ft^2
		g_c = gravitational constant, 32.2 ft/sec^2
		rho = density of propellant, lb/ft^3
		delta_P = pressure drop across orifice, lb/ft^2

	"""
	def propellant_flow_rate(self, A, rho, delta_P, C_d=0.6):
		return (C_d*A)*math.sqrt(2*32.2*rho*delta_P)
	"""

	Design Equation:
		v = C_d * sqrt(2*g_c * (delta_P/rho))

	Variables:
		v = the injection velocity, or velocity of the liquid stream issuing from the orifice, ft/sec
		C_d = orifice discharge coefficient
		g_c = gravitational constant, 32.2 ft/sec^2
		delta_P = pressure drop across orifice, lb/ft^2
		rho = density of propellant, lb/ft^3

	"""
	def injection_velocity(self, delta_P, rho, C_d=0.6):
		return C_d*math.sqrt((2*g_c)*(delta_P/rho))

def main():
    print("Hello World")
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

    w_t = var.total_flow_rate()
    print("total flow rate: " + str(var.total_flow_rate()))
    w_f = var.fuel_flow_rate()
    print("fuel flow rate: " + str(var.fuel_flow_rate()))
    w_o = var.oxidizer_flow_rate()
    print("oxidizer flow rate: " + str(var.oxidizer_flow_rate()))

    noz = nozzle()
    T_t = noz.tempurature_of_gases_at_nozzle_throat(var.T_c, con)
    print("tempurature of gases at nozzle throat: " + str(noz.tempurature_of_gases_at_nozzle_throat(var.T_c, con)))
    P_t = noz.gas_pressure_at_nozzle_throat(var.P_c, con)
    print ("gas presure at nozzle throat: " + str(noz.gas_pressure_at_nozzle_throat(var.P_c, con)))
    A_t = noz.nozzle_throat_cross_sectional_area(w_t, P_t, T_t, con)
    print ("nozzle throat cross-sectional area: " + str(noz.nozzle_throat_cross_sectional_area(w_t, P_t, T_t, con)))
    M_e = noz.mach_number(var.P_c, con)
    print ("mach number: " + str(noz.mach_number(var.P_c, con)))
    A_e = noz.nozzle_exit_cross_sectional_area(A_t, M_e, con)
    print ("nozzle exit cross-sectional area: *broken*" + str(noz.nozzle_exit_cross_sectional_area(A_t, M_e, con)))

    com = combustion_chamber()
    

if __name__ == "__main__":
    main()