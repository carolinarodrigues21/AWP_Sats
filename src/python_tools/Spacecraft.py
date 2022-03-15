'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Spacecraft class definition
'''

# Python standard libraries
import os
import math as m

# 3rd party libraries
import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'dark_background' )
from scipy.integrate import ode
import spiceypy as spice

# AWP libraries
import orbit_calculations as oc
import numerical_tools    as nt
import plotting_tools     as pt
import planetary_data     as pd
import spice_data         as sd

def null_config():
	return {
		'cb'             : pd.earth, #central body
		'date0'          : '2021-04-01', #quer plotar algo que dependa do tempo por ex. posição, rotação, attitude da terra depende do tempo e ai muda a configuração dependendo do tempo
		'et0'            : None, 
		'frame'          : 'J2000', #earth-center propagation standard J2000
		'dt'             : 100, #seconds between orbits states
		'orbit_state'    : [], #fornecer o vetor posição e velocidade forma 1 de cond inicial
		'coes'           : [], #classical orbit elements forma 2 de cond inicial
		'orbit_perts'    : {},
		'propagator'     : 'lsoda', #qual EDO vc quer usar pra propagar essa trajetória
		'stop_conditions': {},
		'print_stop'     : True,
		'mass0'          : 0,
		'output_dir'     : '.',
		'propagate'      : True
	}

class Spacecraft:

	def __init__( self, config ): #config nesse caso é o dicionario de configuração do Spacecraft_hello_world.py
		self.config = null_config() #default parameters do config dicctionary
		#for each key you passed in in sc you overwrite the default
		for key in config.keys():
			self.config[ key ] = config[ key ]

		self.orbit_perts = self.config[ 'orbit_perts' ]
		self.cb          = self.config[ 'cb' ] 

		#se passa coes vai ter que converter para vetor velocidade e posição
		if self.config[ 'coes' ]:
			self.config[ 'orbit_state' ] = oc.coes2state( self.config[ 'coes' ],
				mu = self.config[ 'cb' ][ 'mu' ] )

		#informação pra se for string use como número de períodos 
		if type( self.config[ 'tspan' ] ) == str:
			self.config[ 'tspan' ] = float( self.config[ 'tspan'] ) *\
				oc.state2period( self.config[ 'orbit_state' ], self.cb[ 'mu' ] )

		#numero de steps que a propagação vai fazer em relação ao tspan que passou
		self.steps = int( 
			np.ceil( self.config[ 'tspan' ] / self.config[ 'dt' ] ) + 1 )
		self.step  = 1

		#data arrays, pre-alocate all the memory 
		self.ets    = np.zeros( ( self.steps, 1 ) ) #efemorous times
		self.states = np.zeros( ( self.steps, 7 ) )	#position, velocity and mass
		self.alts   = np.zeros( ( self.steps, 1 ) ) #altitudes

		self.states[ 0, :6 ] = self.config[ 'orbit_state' ]
		self.states[ 0, 6  ] = self.config[ 'mass0'  ]
		self.alts  [ 0 ]     = nt.norm( self.states[ 0, :3 ] ) -\
									self.cb[ 'radius' ]

		self.assign_stop_condition_functions()
		self.assign_orbit_perturbations_functions()
		self.load_spice_kernels()

		if not os.path.exists( self.config[ 'output_dir' ] ):
			os.mkdir( self.config[ 'output_dir' ] )

		self.solver = ode( self.diffy_q )
		self.solver.set_integrator( self.config[ 'propagator' ] )
		self.solver.set_initial_value( self.states[ 0, : ], 0) 

		self.coes_calculated    = False
		self.latlons_calculated = False

		if self.config[ 'propagate' ]:
			self.propagate_orbit()

	def assign_stop_condition_functions( self ):

		self.stop_conditions_map = {
			'max_alt'  : self.check_max_alt,
			'min_alt'  : self.check_min_alt,
			'exit_SOI' : self.check_exit_SOI,
			'enter_SOI': self.check_enter_SOI
			}

		self.stop_condition_functions = [ self.check_deorbit ]

		for key in self.config[ 'stop_conditions' ].keys():
			self.stop_condition_functions.append(
				self.stop_conditions_map[ key ] )

	def assign_orbit_perturbations_functions( self ): #perturbações
	
		self.orbit_perts_funcs_map = {
			'J2': self.calc_J2
		}
		self.orbit_perts_funcs = []

		for key in self.config[ 'orbit_perts' ]:
			self.orbit_perts_funcs.append( 
				self.orbit_perts_funcs_map[ key ] )

	def load_spice_kernels( self ):
		spice.furnsh( sd.leapseconds_kernel )
		self.spice_kernels_loaded = [ sd.leapseconds_kernel ] #spice converte data para seg

		#se passa efemorous time nao-0, em segs,mas senao passar uma data 
		if self.config[ 'et0' ] is not None:
			self.et0 = self.config[ 'et0' ]
		else:
			self.et0 = spice.str2et( self.config[ 'date0' ] )

		self.ets = np.arange( self.et0,
			self.et0 + self.config[ 'tspan' ] + self.config[ 'dt' ],
			self.config[ 'dt' ] )

	def check_deorbit( self ):
		if self.alts[ self.step ] < self.cb[ 'deorbit_altitude' ]:
			if self.config[ 'print_stop' ]:
				self.print_stop_condition( 'deorbit altitude' )
			return False
		return True

	def check_max_alt( self ):
		if self.alts[ self.step ] > self.config[ 'stop_conditions' ][ 'max_alt' ]:
			if self.config[ 'print_stop' ]:
				self.print_stop_condition( 'max altitude' )
			return False
		return True

	def check_min_alt( self ):
		if self.alts[ self.step ] > self.config[ 'stop_conditions' ][ 'min_alt' ]:
			if self.config[ 'print_stop' ]:
				self.print_stop_condition( 'min altitude' )
			return False
		return True

	def check_exit_SOI( self ):
		if nt.norm( self.states[ self.step, :3 ] ) > self.cb[ 'SOI' ]:
			if self.config[ 'print_stop' ]:
				self.print_stop_condition( '%s SOI exit' % self.cb[ 'name' ] )
			return False
		return True

	def check_enter_SOI( self ):
		body      = self.config[ 'stop_conditions' ][ 'enter_SOI' ]
		r_cb2body = spice.spkgps(
						body[ 'SPICE_ID' ], self.ets[ self.step ],
						self.config[ 'frame' ], self.cb[ 'SPICE_ID' ] )[ 0 ]
		r_sc2body = r_cb2body - self.states[ self.step, :3 ]

		if nt.norm( r_sc2body ) < body[ 'SOI' ]:
			self.print_stop_condition( '%s SOI entry' % body[ 'name' ] )
			return False

		return True

	def print_stop_condition( self, parameter ):
		print( f'Spacecraft has reached {parameter}.' )

	def check_stop_conditions( self ):
		for stop_condition in self.stop_condition_functions:
			if not stop_condition():
				return False #se voltar false vai parar a propagação
		return True

	#calcula perturbação 3D em coordenadas cartesianas
	def calc_J2( self, et, state ): 
		z2     = state[ 2 ] ** 2
		norm_r = nt.norm( state[ :3 ] )
		r2     = norm_r ** 2
		tx     = state[ 0 ] / norm_r * ( 5 * z2 / r2 - 1 )
		ty     = state[ 1 ] / norm_r * ( 5 * z2 / r2 - 1 )
		tz     = state[ 2 ] / norm_r * ( 5 * z2 / r2 - 3 )
		return 1.5 * self.cb[ 'J2' ] * self.cb[ 'mu' ] *\
			   self.cb[ 'radius' ] ** 2 \
			 / r2 ** 2 * np.array( [ tx, ty, tz ] )

	# define a ED que mostra como o Spacecraft se move, muda ao longo do tempo
	def diffy_q( self, et, state ): 
		rx, ry, rz, vx, vy, vz, mass = state
		r         = np.array( [ rx,   ry,   rz   ] )
		v         = np.array( [ vx,   vy,   vz   ] )
		norm_r    = nt.norm( r )
		mass_dot  = 0.0 			#derivada da massa no tempo. util no low trust
		state_dot = np.zeros( 7 )	#derivada do estado
		et       += self.et0	

		#calculates the 2-body acceleration
		a = -r * self.cb[ 'mu' ] / norm_r ** 3  

		for pert in self.orbit_perts_funcs:
			a += pert( et, state ) #chama a perturbação J2 aqui

		state_dot[ :3  ] = v # ou [ vx,   vy,   vz   ] 
		state_dot[ 3:6 ] = a
		state_dot[ 6   ] = mass_dot
		return state_dot

	def propagate_orbit( self ):
		print( 'Propagating orbit..' )

		while self.solver.successful() and self.step < self.steps:
			self.solver.integrate( self.solver.t + self.config[ 'dt' ] )
			self.states[ self.step ] = self.solver.y
			self.alts  [ self.step ] = nt.norm( self.solver.y[ :3 ] ) -\
										self.cb[ 'radius' ]
			if self.check_stop_conditions(): #True
				self.step += 1
			else: #False
				break

		self.ets    = self.ets   [ :self.step ]
		self.states = self.states[ :self.step ]
		self.alts   = self.alts  [ :self.step ]

	def calc_coes( self ):
		print( 'Calculating COEs..' )
		self.coes = np.zeros( ( self.step, 6 ) )

		for n in range( self.step ):
			self.coes[ n, : ] = oc.state2coes( 
				self.states[ n, :6 ], mu = self.cb[ 'mu' ] )
			
		self.coes_rel        = self.coes[ : ] - self.coes[ 0, : ]
		self.coes_calculated = True

	def calc_apoapses_periapses( self ): #calculo do apoasis e periasis (afelio e perielio geral)
		if not self.coes_calculated:
			self.calc_coes()

		# outra forma de fazer é -> apoapse = sma * (1 + eccentricity) 
		self.apoapses  = self.coes[ :, 0 ] * ( 1 + self.coes[ :, 1 ] )
		self.periapses = self.coes[ :, 0 ] * ( 1 - self.coes[ :, 1 ] )

		#tudo vetorizado entao vai fornecer um array de apo e peri 
		#util para aerodynamic drag e low trust

	def plot_3d( self, args = { 'show': True } ):
		pt.plot_orbits( [ self.states[ :, :3 ] ], args )

	def plot_groundtracks( self, args = { 'show': True } ):
		if not self.latlons_calculated:
			self.calc_latlons()

		pt.plot_groundtracks( [ self.latlons[ : ] ], args )

	def plot_coes( self, args = { 'show': True }, step = 1 ):
		if not self.coes_calculated:
			self.calc_coes()

		pt.plot_coes( self.ets[ ::step ], self.coes[ ::step ], args )
