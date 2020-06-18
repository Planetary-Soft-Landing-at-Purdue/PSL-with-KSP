from multiprocessing import Process, Value, Array
from config import *
import krpc, time, pdg
import numpy as np


def process_time():
	global conn, sc, vessel, orbit, body, bcbf, pcpf
	global position_stream, velocity_stream, mass_stream, met_stream

	conn = krpc.connect(address='192.168.1.181')
	sc = conn.space_center
	vessel = sc.active_vessel
	orbit = vessel.orbit
	body = orbit.body
	bcbf = body.reference_frame
	pcpf = conn.space_center.ReferenceFrame.create_relative(
		bcbf,
		position=vessel.position(bcbf),
		rotation=vessel.rotation(bcbf))

	position_stream = conn.add_stream(vessel.position, pcpf)
	position_stream.add_callback(position_callback)
	position_stream.start()

	velocity_stream = conn.add_stream(vessel.velocity, pcpf)
	velocity_stream.add_callback(velocity_callback)
	velocity_stream.start()
	
	mass_stream = conn.add_stream(getattr, vessel, 'mass')
	mass_stream.add_callback(mass_callback)
	mass_stream.start()
	
	met_stream = conn.add_stream(getattr, vessel, 'met')
	met_stream.add_callback(met_callback)
	met_stream.start()

	while True: pass


def met_callback(self):
	met.value = self
	control()
def mass_callback(self): mass.value = self
def position_callback(self):
	position[0].value = self[0]
	position[1].value = self[1]
	position[2].value = self[2]
def velocity_callback(self):
	velocity[0].value = self[0]
	velocity[1].value = self[1]
	velocity[2].value = self[2]

def control():
	global cet0
	if new_eta.value == 1:
		cet0 = met.value
		cet = 0
		new_eta.value = 0
	cet = met.value-cet0
	n = 4*int(cet/pdg.delta_t)
	vessel.auto_pilot.direction = eta[n].value, eta[n+1].value, eta[n+2].value
	vessel.control.throttle = eta[n+3].value


def process_guid(): pass


if __name__ == '__main__':

	process_time = Process(target=process_time,args=(),name='TIME')
	# process_guid = Process(target=process_guid,args=(),name='GUID')

	process_time.start()
	# process_guid.start()




	global tWait, tSolve, tDist, eta
	while met.value == 0: pass

	state0 = np.array([
		position[0].value, position[1].value, position[2].value,
		velocity[0].value, velocity[1].value, velocity[2].value,
		np.log(mass.value), 0, 0, 0, 0])
	tWait, tSolve, tDist = pdg.findPath(pdg.delta_t, state0, initialSearch=True)

	while tSolve > 0:
		t0 = met.value
		path = pdg.findPath(pdg.delta_t, state0, tWait=tWait, tSolve=tSolve, tDist=tDist)
		eta = Array('d', path)
		tSolve -= met.value-t0
		new_eta.value = 1

	process_time.terminate()