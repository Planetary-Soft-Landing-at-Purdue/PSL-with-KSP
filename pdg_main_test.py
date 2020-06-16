from multiprocessing import Process, Manager
from multiprocessing.managers import Namespace
from time import sleep, time
import krpc


#  Connects to KRPC server and updates Namespace
#	of KSP attributes on each KSP timestep.
def process_time():
	global conn, vessel, orbit, body, bcbf
	conn = krpc.connect(address='192.168.1.169')
	vessel = conn.space_center.active_vessel
	orbit = vessel.orbit
	body = orbit.body
	bcbf = body.reference_frame

	met_stream = conn.add_stream(getattr, vessel, 'met')
	met_stream.add_callback(timestep)
	met_stream.start()

	while True: pass


#  Performs all actions needed on each KSP timestep.
def timestep(self):

	#  Update state.
	ns.position = vessel.position(bcbf)
	ns.velocity = vessel.velocity(bcbf)
	ns.mass = vessel.mass
	t0 = time()
	ns.met = self
	print(time()-t0)

	#  Run control function.
	# control()


#  Sends control commands to vessel.
def control():
	global cet0
	while not hasattr(ns, 'new_eta'): pass
	if ns.new_eta == True:
		cet0 = ns.met
		cet = 0
		ns.new_eta = False
	cet = ns.met-cet0
	print('CET:', round(cet, 3), 'CET0:', round(cet0, 3))


def process_guid():
	while not hasattr(ns, 'met'): pass
	while True:
		ns.new_eta = True
		sleep(5)


if __name__ == '__main__':

	manager = Manager()
	ns = manager.Namespace()

	Process_time = Process(target=process_time)
	Process_guid = Process(target=process_guid)
	Process_time.start()
	Process_guid.start()
	Process_time.join()
	Process_guid.join()