'''
Print Symbolic Euler Angle Rotation Matrix


Calcular a matriz de cossenos de direção associada
com a sequencia 313 dos angulos de euler. Referente
a um eixo inercial. 
'''

#importanto as bibliotecas usadas 
from sys import path 
path.append ('../AWP/src/python_tools')
import plotting_tools as pt
import spiceypy as spice
import math as mt


#conversão para radianos
d2r = mt.pi / 180.0 

#nome para a imagem
fn  = 'frames.png'

if __name__ == '__main__':

	#eixo inercial 
	eul_iniercial = spice.eul2m( 0, 0, 0, 3, 1, 1) 

	#angulos que caracterizam uma orbita Kepleriana
	raan   =  0  * d2r	#right ascension
	inc    =  90  * d2r #inclinacao
	aop    =  60 * d2r	#argumento de periapsis

	eul313 = spice.eul2m( aop, inc, raan, 3, 1, 1 )
	#frames = [ eul313, eul313.T ]
	frames = [ eul313]

	print( 'eul2m' )
	print( eul313 )
	print()
	print( 'transpose' )
	print( eul313.T )

	config = {
		'frame_colors': [ 'm', 'c' ],
		#'frame_labels': [ 'eul2m', 'Transpose' ],
		'frame_labels': [ 'eul2m'],
		'elevation'   : 12,
		'azimuth'     : -19,
		#'filename'    : fn
		'show': True
	}

	pt.plot_reference_frames( frames, config )