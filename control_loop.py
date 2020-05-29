import krpc
import numpy as np

def inner_loop(yeet):
	global pdg_t0							#  refer to global variable
	curr_met = met()						#  get current MET
	pdg_et = curr_met - pdg_t0				#  update PDG elapsed time
	p_index = p[int(round(pdg_et/pdg_dt))]	#  set p matrix index
	#  if met has grown x seconds, run outer loop

def outer_loop():
	pass	#  run pdg, update p matrix

#  connect and set active vessel, IPv4
conn = krpc.connect(address = '192.168.1.169')
sc = conn.space_center
vessel = sc.active_vessel
ap = vessel.auto_pilot

#  timestep and p matrix
pdg_dt = 1
p = np.zeros((   69   , 4))		#  69 placeholder

#  create met stream
met = conn.add_stream(getattr, vessel, 'met')
print('MET stream done')

#  create and start MET callback
pdg_t0 = met()		#  set PDG t0
met_last = pdg_t0	#  set last MET for outer loop
met.add_callback(inner_loop)
met.start()

while True:
	pass