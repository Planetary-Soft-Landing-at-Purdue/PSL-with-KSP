from multiprocessing import Process, Manager, Lock
import krpc, time, pdg, math
import numpy as np

def process_time(ns):
	global conn, sc, vessel, orbit, body, bcbf, bci, omega, pcpf

	ns.startPDG = False

	# makes all necessary connections with ksp
	conn 	= krpc.connect(address='128.211.208.159')
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
	met_stream 		= conn.add_stream(getattr, vessel, 'met')
	position_stream.start()
	velocity_stream.start()
	mass_stream.start()
	met_stream.start()

	# launches the vessel, waits until it reaches a certain 
	# altitude, then calls pdg and begins descent
	print("Made connections, launching")

	vessel.auto_pilot.sas = True
	vessel.control.throttle = .6
	vessel.control.activate_next_stage()

	while position_stream()[1] < 400: pass

	vessel.control.throttle = 0
	vessel.auto_pilot.sas 	= False
	vessel.auto_pilot.engage()
	vessel.auto_pilot.reference_frame = pcpf

	print("Starting controlled descent")
	startTime, n, tPrev, ns.eta = time.time(), 0, -ns.dt, [0, 0, 0, 0]

	# continuously updates the vessel's thrust until it reaches a certain
	# altitude above its target landing spot
	while ns.new_eta == None or position_stream()[1] > 1:
		# updates current state vector
		position, velocity, mass = position_stream(), velocity_stream(), mass_stream()
		ns.stateVect = np.array([
			 position[1], position[0], position[2],
			 velocity[1], velocity[0], velocity[2],
			 np.log(mass), ns.eta[n]*ns.eta[n+3], ns.eta[n+1]*ns.eta[n+3], ns.eta[n+2]*ns.eta[n+3], ns.eta[n+3]])
		ns.startPDG = True

		# waits until a new eta has been found
		if ns.new_eta == None: 
			vessel.control.throttle = .01
			time.sleep(ns.dt)

		else:
			if ns.new_eta == True:
				ns.new_eta = False
				startTime, n, tPrev = time.time(), 0, -ns.dt

			# waits until a there's a new discrete time step before updating the
			# vessel's thrust
			while int((time.time() - startTime) / ns.dt) == tPrev: pass
			n, tPrev = n + 4, int((time.time() - startTime) / ns.dt)

			vessel.auto_pilot.target_direction = ns.eta[n+1], ns.eta[n], ns.eta[n+2]
			vessel.control.throttle 		   = ns.eta[n+3] * (mass_stream() / pdg.rho_2)

	ns.startPDG = False
	vessel.control.throttle = 0

def process_guid(ns):
	global tWait, tSolve, dMax, dt

	ns.dt 		 = 1
	ns.new_eta 	 = None

	time.sleep(1)
	while ns.startPDG == False: pass

	timeWait = pdg.PDG(ns.dt, ns.stateVect, initialSearch=True)
	print("Will wait", timeWait, "seconds before starting pdg")
	time.sleep(timeWait)
	print("Stopped waiting")

	ns.dt 		    = .2
	tSolve, dMax, _ = pdg.PDG(ns.dt, ns.stateVect, minDistance=True)
	tSolveTotal     = tSolve

	while ns.startPDG and tSolve > 0:
		startTime = time.time()

		ns.eta, _  = pdg.PDG(ns.dt, ns.stateVect, tSolveTotal=tSolveTotal, tWait=timeWait, tSolve=tSolve, dMax=dMax)
		print("Found new solution", tSolve)
		ns.new_eta = True

		time.sleep(1)
		tSolve -= time.time() - startTime
	
		if tSolve < tSolveTotal / 4 : ns.dt = .1
		#if tSolve < tSolveTotal / 8 : ns.dt = .05
		#if tSolve < tSolveTotal / 8: ns.dt = .02

	#process_time.terminate()

if __name__ == '__main__':
	lock = Lock()
	manager = Manager()
	ns = manager.Namespace()

	#initiating and starting two processes
	Process_time = Process(target=process_time, name='TIME', args=(ns,))
	Process_guid = Process(target=process_guid, name='GUID', args=(ns,))

	Process_time.start()
	Process_guid.start()
	Process_time.join()
	Process_guid.join()


