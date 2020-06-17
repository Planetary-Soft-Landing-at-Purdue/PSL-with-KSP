from multiprocessing import Process, Value, Lock
from config import *
import krpc, time, pdg
import numpy as np

def process_time():
	global conn, vessel, orbit, body, bcbf
	global position_stream, velocity_stream, mass_stream, met_stream
	conn = krpc.connect(address='192.168.1.181')
	vessel = conn.space_center.active_vessel

	#  Create reference frame at launch pad
	pcpf = 

	position_stream = conn.add_stream(vessel.position, bcbf)
	velocity_stream = conn.add_stream(vessel.velocity, bcbf)
	mass_stream = conn.add_stream(getattr, vessel, 'mass')
	met_stream = conn.add_stream(getattr, vessel, 'met')

	met_stream.add_callback(timestep)

	position_stream.start()
	velocity_stream.start()
	mass_stream.start()
	met_stream.start()
	while True: pass


def timestep(self):
	met.value = self
	update_state()
	control()


def update_state():
	with lock:
		position = position_stream()
		velocity = velocity_stream()
		mass = mass_stream()



def control():
	global cet0
	if new_eta.value == 1:
		cet0 = met.value
		cet = 0
		new_eta.value = 0
	cet = met.value-cet0
	with lock: print('  CET:', round(cet, 3), 'MET:', round(met.value, 3))
	

def process_guid():
	
	#  Wait for time counter to begin
	while met.value == 0: pass

	#  Set initial state
	state0 = np.array([
		position[0], position[1], position[2],
		velocity[0], velocity[1], velocity[2],
		np.log(mass), 0, 0, 0, 0])
	tWait, tSolve, tDist = pdg.findPath(pdg.delta_t, state0,)
	while solvetime > 0:
		t0 = met.value
		time.sleep(1)
		d_solvetime = met.value-t0
		with lock:
			print('GdT:', round(d_solvetime, 3), ';', 'solvetime:', round(solvetime, 2))
		solvetime -= d_solvetime
		new_eta.value = 1
	process_time.terminate()


if __name__ == '__main__':

	lock = Lock()

	process_time = Process(target=process_time,args=(),name='TIME')
	process_guid = Process(target=process_guid,args=(),name='GUID')

	process_time.start()
	process_guid.start()