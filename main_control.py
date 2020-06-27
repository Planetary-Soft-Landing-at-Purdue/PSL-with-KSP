from multiprocessing import Process, Manager, Lock
import krpc, time, pdg, math
import numpy as np

def process_time():
	global conn, sc, vessel, orbit, body, bcbf, bci, omega, pcpf

	ns.startPDG = False
	targetHeight = 20

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

	vessel.control.throttle = 0
	vessel.auto_pilot.sas 	= False
	vessel.auto_pilot.engage()
	vessel.auto_pilot.reference_frame = pcpf

	print("Starting controlled descent")
	startTime, n, tPrev, ns.eta = time.time(), 0, -ns.dt, [0, 0, 0, 0]

	# continuously updates the vessel's thrust until it reaches a certain
	# altitude above its target landing spot
	while ns.new_eta == None or position_stream()[1] > targetHeight + 5:
		# waits until a new eta has been found
		position, velocity, mass = position_stream(), velocity_stream(), mass_stream()
		ns.stateVect = np.array([
			 position[1] - targetHeight, position[0], position[2],
			 velocity[1]		       , velocity[0], velocity[2],
			 np.log(mass), ns.eta[n]*ns.eta[n+3], ns.eta[n+1]*ns.eta[n+3], ns.eta[n+2]*ns.eta[n+3], ns.eta[n+3]])
		ns.startPDG = True

		while ns.new_eta == None: pass 
		# updates current state vector
		
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
	
	# makes thrust in opposite direction of horizontal velo until horizontal velo reaches 0
	print("Stabilizing horizontal velocities")
	prop = .15
	while position_stream()[1] > 5 and (abs(velocity_stream()[0]) > .1 or abs(velocity_stream()[2]) > .1):
		dx, dz, dy 			= velocity_stream()[0], velocity_stream()[1], velocity_stream()[2]
		dir_x, dir_z, dir_y	= -prop * dx, -dz, -prop * dy

		if (dir_x**2 + dir_y**2)**.5 / dir_z > math.tan(pdg.theta):
			dir_z = (dir_x**2 + dir_y**2)**.5 / math.tan(pdg.theta)

		dir_mag	= (dir_x**2 + dir_z**2 + dir_y**2) ** .5 

		vessel.auto_pilot.target_direction 	= (dir_x / dir_mag), (dir_z / dir_mag), (dir_y / dir_mag)
		vessel.control.throttle 			= 9.81 * (mass_stream() / pdg.rho_2)	

	# makes a positive thrust until the vessel's vertical velocity is zero
	print("Starting vertical descent")
	vessel.auto_pilot.target_direction = 0, 1, 0
	while velocity_stream()[1] < -1.0: 
		vessel.control.throttle = 2 * 9.81 * (mass_stream() / pdg.rho_2)		
	# has a positive thrust slightly less than gravity until it lands
	vessel.control.throttle = .7 * 9.81 * (mass_stream() / pdg.rho_2)		
	time.sleep(1)
	while position_stream()[1] > 2:
		vessel.control.throttle = .95 * 9.81 * (mass_stream() / pdg.rho_2)		
	
	vessel.control.throttle = 0

def process_guid():
	global tWait, tSolve, dMax, dt

	ns.dt 		 = .5
	ns.new_eta 	 = None

	while ns.startPDG == False: pass

	while ns.startPDG:
		tWait, tSolve, dMax = pdg.PDG(ns.dt, ns.stateVect, initialSearch=True)
		if int(tSolve / ns.dt) > 5:
			ns.eta, _   = pdg.PDG(ns.dt, ns.stateVect, tWait=tWait, tSolve=tSolve, dMax=dMax)
			print("Found solution")
			ns.new_eta 	= True
		time.sleep(5)

	#process_time.terminate()

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


