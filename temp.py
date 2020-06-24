from multiprocessing import Process, Manager, Lock
import krpc, time, pdg
import numpy as np

def process_time():
	global conn, sc, vessel, orbit, body, bcbf, bci, omega, pcpf

	ns.startPDG = False

	# makes all necessary connections with ksp
	conn 	= krpc.connect(address='192.168.0.109')
	sc 		= conn.space_center
	vessel 	= sc.active_vessel
	orbit 	= vessel.orbit
	body 	= orbit.body
	bcbf 	= body.reference_frame
	bci 	= body.non_rotating_reference_frame
	omega 	= body.angular_velocity(bci)

	# creates reference frame, with origin at the launching pad
	pcpf = sc.ReferenceFrame.create_relative( bcbf,
		position=vessel.position(bcbf),
		rotation=vessel.rotation(bcbf))

	# streams are set up to create efficient flow of data from ksp
	position_stream = conn.add_stream(vessel.position, pcpf)
	velocity_stream = conn.add_stream(vessel.velocity, pcpf)
	mass_stream 	= conn.add_stream(getattr, vessel, 'mass')
	position_stream.start()
	velocity_stream.start()
	mass_stream.start()

	# launches the vessel, waits until it reaches a certain 
	# altitude, then calls pdg and begins descent
	print("Made connections, launching")

	vessel.auto_pilot.sas = True
	vessel.control.throttle = .6
	vessel.control.activate_next_stage()

	while position_stream()[1] < 800: pass

	position, velocity, mass = position_stream(), velocity_stream(), mass_stream()
	ns.stateVect = np.array([
			 position[1], position[0], position[2],
			-velocity[1], velocity[0], velocity[2],
			 np.log(mass), 0, 0, 0, 0])

	ns.startPDG = True
	vessel.control.throttle = 0
	vessel.auto_pilot.sas 	= False
	vessel.auto_pilot.engage()
	vessel.auto_pilot.reference_frame = pcpf

	# waits until the vessel comes back to the same altitude
	# at which it cut off its engine
	while position_stream()[1] > 800: pass

	print("Starting controlled descent")
	startTime, n, tPrev = time.time(), 0, -ns.dt

	while position_stream()[1] > 50:
		# waits until a new eta has been found
		while ns.new_eta == None: pass 
		if ns.new_eta == True:
			ns.new_eta = False
			startTime, n, tPrev = time.time(), 0, -ns.dt
		# waits until a there's a new discrete time step before updating the
		# vessel's thrust
		while int((time.time() - startTime) / ns.dt) == tPrev: pass
		n, tPrev = n + 4, int((time.time() - startTime) / ns.dt)

		vessel.auto_pilot.target_direction = ns.eta[n+1], ns.eta[n], ns.eta[n+2]
		vessel.control.throttle 		   = ns.eta[n+3] * (mass_stream() / pdg.rho_2)

		# updates current state vector
		position, velocity, mass = position_stream(), velocity_stream(), mass_stream()
		ns.stateVect = np.array([
			 position[1], position[0], position[2],
			 velocity[1], velocity[0], velocity[2],
			 np.log(mass), ns.eta[n]*ns.eta[n+3], ns.eta[n+1]*ns.eta[n+3], ns.eta[n+2]*ns.eta[n+3], ns.eta[n+3]])

	print("Starting vertical descent")

	vessel.auto_pilot.target_direction = 0, 1, 0
	vessel.control.throttle 		   = .5 * 9.81 * (mass_stream() / pdg.rho_2)
	time.sleep(1)

	while position_stream()[1] > 5:
		vessel.auto_pilot.target_direction = 0, 1, 0
		vessel.control.throttle 		   = .95 * 9.81 * (mass_stream() / pdg.rho_2)

	print("Cutting engines")
	vessel.control.throttle = 0

def process_guid():
	global tWait, tSolve, dMax, dt

	ns.dt 		 = .75
	ns.new_eta 	 = None

	while ns.startPDG == False: pass

	while True:
		tWait, tSolve, dMax = pdg.PDG(ns.dt, ns.stateVect, initialSearch=True)
		ns.eta, _   = pdg.PDG(ns.dt, ns.stateVect, tWait=tWait, tSolve=tSolve, dMax=dMax)
		ns.new_eta 	= True
		print("Found new solution")
		time.sleep(5)

	process_time.terminate()

if __name__ == '__main__':
	lock = Lock()
	manager = Manager()
	ns = manager.Namespace()

	#initiating and starting two processes
	Process_time = Process(target=process_time, name='TIME')
	Process_guid = Process(target=process_guid, name='GUID')

	Process_time.start()
	Process_guid.start()
	Process_time.join()
	Process_guid.join()


