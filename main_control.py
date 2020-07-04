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

	# launches the vessel, waits until it reaches a certain 
	# altitude, then calls pdg and begins descent
	print("Made connections, launching")

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

	realSol, pdgSol = [], []
	ns.predTime     = 1

	# continuously updates the vessel's thrust until it reaches a certain
	# altitude above its target landing spot
	while ns.new_eta == None or position_stream()[1] > 1:
		# updates ns namespace variables
		pos, velo, mass = position_stream(), velocity_stream(), mass_stream()
		ns.predPos  = [pos[0]  + ns.predTime * velo[0],     pos[1]  + ns.predTime * velo[1],               pos[2]  + ns.predTime * velo[2]     ]
		ns.predVelo = [velo[0] + ns.predTime * ns.eta[n+1], velo[1] + ns.predTime * (ns.eta[n] + ns.g[7]), velo[2] + ns.predTime * ns.eta[n+2] ]
		ns.nextEta  = ns.eta[n], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3] 
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

				pdgSol.extend(ns.sol[(n // 4) * 11 : (n // 4) * 11 + 11])
				realSol.extend([position_stream()[1], position_stream()[0], position_stream()[2], 
					velocity_stream()[1], velocity_stream()[0], velocity_stream()[2], np.log(mass_stream()),
					ns.eta[n], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3]])

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

	while ns.startPDG:
		startTime = ns.met

		# update state vector
		pdg.update_state_vector(ns.predPos, ns.predVelo, ns.predMass, ns.nextEta)

		# calls pdg to optimize fuel use and recalculate path, tells process time
		# that there is a new solution
		ns.sol, ns.eta = pdg.PDG(tSolve=tSolve)
		print("----", tSolve)

		ns.new_eta  = True
		tSolve     -= ns.met - startTime
		ns.predTime = ns.met - startTime
	
		#if tSolve < tSolveTotal / 4 : ns.dt = .1
		#if tSolve < tSolveTotal / 8 : ns.dt = .05
		#if tSolve < tSolveTotal / 8: ns.dt = .02

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