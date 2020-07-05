from multiprocessing import Process, Manager, Lock
from math import log, pi, sin, cos
from pdg import PDG
import krpc, time, math
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
	omega 	= -1*body.angular_velocity(bci)[1]

	# creates reference frame, with origin at the launching pad
	pcpf = sc.ReferenceFrame.create_relative( bcbf,
		position=vessel.position(bcbf),
		rotation=vessel.rotation(bcbf))

	# defines rotational velocity and the gravitation acceleration of the planet
	lat = body.latitude_at_position(vessel.position(bcbf), bcbf) * pi / 180
	lon = body.longitude_at_position(vessel.position(bcbf), bcbf) * pi / 180
	ns.w = [omega * sin(lat + pi / 2), omega * cos(lat + pi / 2) * cos(lon), omega * cos(lat + pi / 2) * sin(lon)]
	ns.g = [0, 0, 0, 0, 0, 0, 0, -body.surface_gravity, 0, 0, 0]
	
	# streams are set up to create efficient flow of data from ksp
	position_stream = conn.add_stream(vessel.position, pcpf)
	velocity_stream = conn.add_stream(vessel.velocity, pcpf)
	mass_stream 	= conn.add_stream(getattr, vessel, 'mass')
	met_stream 		= conn.add_stream(getattr, vessel, 'met')
	position_stream.start()
	velocity_stream.start()
	mass_stream.start()
	met_stream.start()

	print("Made connections, launching")

	# launches the vessel, waits until it reaches a certain altitude, then cuts engines 
	# and waits until vertical velocity hits zero, then starts guided descent
	vessel.auto_pilot.sas = True
	vessel.control.throttle = .6
	vessel.control.activate_next_stage()

	while position_stream()[1] < 400: pass
	vessel.control.throttle = 0
	while velocity_stream()[1] > 0: pass

	vessel.auto_pilot.sas 	= False
	vessel.auto_pilot.engage()
	vessel.auto_pilot.reference_frame = pcpf

	print("Starting controlled descent")
	startTime, n, tPrev, ns.eta = met_stream(), 0, -1, [0, 0, 0, 0]

	# realSol and pdgSol are lists that keep track of the actual state vector at every time step and pdg-predicted 
	# state vector at every time step. predTime is the time offset from the when pdg is called and when it finishes solving.
	realSol, pdgSol    = [], []
	ns.predTime, ns.dt = 1, 1

	# continuously updates the vessel's thrust until it reaches a certain
	# altitude above its target landing spot
	while ns.new_eta == None or position_stream()[1] > 2:
		# predicts what the state vector will be when the next pdg call finishes, pdg solves for this state vector
		pos, velo, mass = position_stream(), velocity_stream(), mass_stream()
		ns.predPos  = [pos[0]  + ns.predTime * velo[0],     pos[1]  + ns.predTime * velo[1],               pos[2]  + ns.predTime * velo[2]     ]
		ns.predVelo = [velo[0] + ns.predTime * ns.eta[n+1], velo[1] + ns.predTime * (ns.eta[n] + ns.g[7]), velo[2] + ns.predTime * ns.eta[n+2] ]
		#ns.nextEta  = ns.eta[n], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3] 
		# if there is a next eta, send the next eta, if no next eta, send current eta
		nextEta = n + 4 * int(ns.predTime / ns.dt)
		if nextEta < len(ns.eta):
			ns.nextEta = ns.eta[nextEta: nextEta + 4]
		else:
			ns.nextEta = ns.eta[-4:]

		ns.predMass = mass_stream()
		ns.met      = met_stream()

		ns.startPDG = True

		# waits until a new eta has been found
		if ns.new_eta == None: pass
		else:
			# if Guidance found a new path, reset startTime, n, and tPrev
			if ns.new_eta == True:
				ns.new_eta = False
				startTime, n, tPrev = met_stream(), 0, -1

			# if an entire time step has passed, move to next set of eta values
			if int((met_stream() - startTime) / ns.dt) != tPrev: 
				n, tPrev = n + 4, int((met_stream() - startTime) / ns.dt)

				# adds real state vector and pdg state vector to realSol and pdgSol
				pdgSol.extend(ns.sol[(n // 4) * 11 : (n // 4) * 11 + 11])
				realSol.extend([position_stream()[1], position_stream()[0], position_stream()[2], 
					velocity_stream()[1], velocity_stream()[0], velocity_stream()[2], np.log(mass_stream()),
					ns.eta[n], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3]])

			# sets vessel direction and throttle based on current eta
			vessel.auto_pilot.target_direction = ns.eta[n+1], ns.eta[n], ns.eta[n+2]
			vessel.control.throttle 		   = ns.eta[n+3] * (mass_stream() / ns.rho_2)

	print("Finished Guidance")
	ns.startPDG = False
	vessel.control.throttle = 0

	realData = open("real_data.txt", "w")
	realData.write(str(realSol))
	realData.close()

	pdgData = open("pdg_data.txt", "w")
	pdgData.write(str(pdgSol))
	pdgData.close()

def process_guid():
	global tWait, tSolve, dMax, dt

	time.sleep(1)
	ns.new_eta = None

	# vessel specific constants
	m_f   = 9495
	rho_1 = 936508 * .15
	rho_2 = 936508 * .4
	alpha = 3.46e-4
	gamma = pi / 3
	theta = pi / 4

	pdg = PDG(m_f, rho_1, rho_2, alpha, gamma, theta)

	# waits until process time tells it to start looking for a solution
	while ns.startPDG == False: pass

	pdg.init_constants(ns.g, ns.w)
	pdg.dt, ns.dt, ns.rho_2 = 1, 1, pdg.rho_2
	pdg.update_state_vector(ns.predPos, ns.predVelo, ns.predMass, ns.nextEta)

	tSolve        = pdg.PDG(findTimeSolve=True)
	print("Found tSolve")
	pdg.dt, ns.dt = .1, .1
	initialSol    = True

	while ns.startPDG:
		# calls pdg to optimize fuel use and recalculate path, tells process time
		# that there is a new solution

		startTime = ns.met

		# update state vector
		pdg.update_state_vector(ns.predPos, ns.predVelo, ns.predMass, ns.nextEta)

		# if time runs out, add two seconds to tSolve
		if tSolve < 1.0: tSolve += 2.0
		# if tSolve is less than 6.0, replace min distance constraint with 2x current
		# distance, allows the vessel to do a vertical descent landing
		if tSolve < 6.0:
			pdg.dMax = 2 * np.linalg.norm(pdg.x[1:3])

		ns.sol, ns.eta = pdg.PDG(tSolve=tSolve)
		print("----", tSolve, pdg.dMax)

		# tells process time that there is a new solution, decrement tSolve based on how much 
		# time has passed, and adjust predTime based on how long the prev pdg took to solve
		tSolve     -= ns.met - startTime
		ns.predTime = ns.met - startTime
		ns.new_eta  = True

		if initialSol:
			initPdg = ns.sol
			initialSol = False

	initPdgData = open("initial_pdg_data.txt", "w")
	pdgStr = ''
	for s in initPdg:
		pdgStr += str(s) + ","
	initPdgData.write(pdgStr)
	initPdgData.close()

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