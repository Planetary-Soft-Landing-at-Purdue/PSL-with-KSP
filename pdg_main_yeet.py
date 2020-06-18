from multiprocessing import Process, Manager
import krpc, time, pdg
import numpy as np


def process_time():
	conn = krpc.connect(address='192.168.1.181')
	sc = conn.space_center
	vessel = sc.active_vessel
	orbit = vessel.orbit
	body = orbit.body
	bcbf = body.reference_frame
	bci = body.non_rotating_reference_frame
	omega = body.angular_velocity(bci)

	pcpf = sc.ReferenceFrame.create_relative(
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
	ns.met = self
	control()
def mass_callback(self): ns.mass = self
def position_callback(self): ns.position = self
def velocity_callback(self): ns.velocity = self

def control():
	global cet0
	if ns.new_eta == True:
		cet0 = ns.met
		cet = 0
		ns.new_eta = False
	cet = ns.met-cet0

	n = 4 * int(cet/pdg.delta_t)
	vessel.auto_pilot.direction = ns.eta[n+2], ns.eta[n], ns.eta[n+1]
	vessel.control.throttle = ns.eta[n+3]

def process_guid():
	global tWait, tSolve, tDist, dt
	dt = 0.5

	while ns.met == 0: pass
	tWait, tSolve, dMax = pdg.findPath(dt, state(ns), initialSearch=True)

	while tSolve > 0:
		t0 = ns.met
		x0 = state(ns)
		ns.eta = pdg.findPath(dt, x0, tWait=tWait, tSolve=tSolve, tDist=tDist)
		tSolve -= met.value-t0
		ns.new_eta = True
	process_time.terminate()

def state(ns):
	return np.array([
		ns.position[0], ns.position[1], ns.position[2],
		ns.velocity[0], ns.velocity[1], ns.velocity[2],
		np.log(ns.mass), 0, 0, 0])


if __name__ == '__main__':
	manager = Manager()
	ns = manager.Namespace()
	Process_time = Process(target=process_time, name='TIME')
	Process_guid = Process(target=process_guid, name='GUID')
	Process_time.start()
	Process_guid.start()
	Process_time.join()
	Process_guid.join()