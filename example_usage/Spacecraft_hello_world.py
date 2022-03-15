'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Hello world of Spacecraft class
Two-body propagation with J2 perturbation for 100 periods
'''

# Python standard libraries
from sys import path 
path.append ('../src/python_tools')

# AWP libraries
from Spacecraft import Spacecraft as SC
from planetary_data import earth

if __name__ == '__main__':
	# condicoes iniciais da Spacecraft: [semi-major axis, eccentricity, inclination in degrees,]
	coes = [ earth[ 'radius' ] + 1000, 0.05, 30.0, 0.0, 0.0, 0.0 ]	#elementos orbitais de kepler
	sc   = SC(
			{
			'coes'       : coes,	#initial conditions of the spacecraft
			'tspan'      : '100', 	#'100' fica com o band	#how long do you want to propagate this orbit for - string: numero de períodos
			'dt'         : 100.0,	#time steps in which you want to see the orbital states are
			'orbit_perts': { 'J2': True }, #perturbações 
			#'propagator' : 'dopris',
			#'propagate'  : False
			} )
	#acima é o dicionario de confuguração pra dizer como vc quer ver a spacecraft 
	sc.plot_3d() 