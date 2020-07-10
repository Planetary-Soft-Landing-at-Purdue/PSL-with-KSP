from multiprocessing import Process, Manager, Lock
from math import log, pi, sin, cos
from pdg import PDG
import krpc, time, math
import numpy as np

def process_vess():	
	ns.startPDG = False

	# makes all necessary connections with ksp
	conn 	= krpc.connect(address='192.168.0.109')
	sc 		= conn.space_center
	vessel 	= sc.active_vessel
	orbit 	= vessel.orbit
	body 	= orbit.body
	bcbf 	= body.reference_frame
	bci 	= body.non_rotating_reference_frame
	omega 	= -1 * body.angular_velocity(bci)[1]

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

	ns.predTime, ns.thrustVector, ns.nextEta = 1, [0, 0, 0, 0], [0, 0, 0, 0]

	while ns.startPDG == False or ns.pos[1] > 2:
		ns.pos, ns.velo, ns.mass = position_stream(), velocity_stream(), mass_stream()
		ns.predPos  = [ns.pos[0]  + ns.predTime * ns.velo[0], 
					   ns.pos[1]  + ns.predTime * ns.velo[1],            
					   ns.pos[2]  + ns.predTime * ns.velo[2]]
		ns.predVelo = [ns.velo[0] + ns.predTime *  ns.thrustVector[1]           ,  
					   ns.velo[1] + ns.predTime * (ns.thrustVector[0] + ns.g[7]), 
					   ns.velo[2] + ns.predTime *  ns.thrustVector[2]           ]
 
		ns.predMass = mass_stream()
		ns.met      = met_stream()

		ns.startPDG = True

		vessel.auto_pilot.target_direction = ns.thrustVector[1], ns.thrustVector[0], ns.thrustVector[2]
		vessel.control.throttle 		   = ns.thrustVector[3] * (mass_stream() / ns.rho_2)

	vessel.control.throttle = 0
	ns.startPDG = False

def process_time():
	time.sleep(1)
	# realSol and pdgSol are lists that keep track of the actual state vector at every time step and pdg-predicted 
	# state vector at every time step. predTime is the time offset from the when pdg is called and when it finishes solving.
	while ns.startPDG == False or ns.new_eta == None: pass

	realSol, pdgSol = [], []

	while ns.pos[1] > 2:
		# if Guidance found a new path, reset t_start, n, and tPrev
		if ns.new_eta == True:
			ns.new_eta = False
			t_start    = ns.met

		n = 4 * int((ns.met - t_start) / ns.dt)

		# sets vessel direction and throttle based on current eta
		ns.thrustVector = ns.eta[n : n + 4]
		if n + 7 < len(ns.eta):
			ns.nextEta  = ns.eta[n + 4 : n + 8]
		else:
			ns.nextEta  = ns.eta[n : n + 4]

		# adds real state vector and pdg state vector to realSol and pdgSol
		pdgSol.extend(ns.sol[(n // 4) * 11 : (n // 4) * 11 + 11])
		realSol.extend([ns.pos[1], ns.pos[0], ns.pos[2], 
			ns.velo[1], ns.velo[0], ns.velo[2], np.log(ns.mass),
			ns.eta[n], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3]])

	realData = open("data/real_data.txt", "w")
	realData.write(str(realSol))
	realData.close()

	pdgData = open("data/pdg_data.txt", "w")
	pdgData.write(str(pdgSol))
	pdgData.close()

def process_guid():
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

	pdg.INIT_CONSTANTS(ns.g, ns.w)
	pdg.dt, ns.dt, ns.rho_2 = 1, 1, pdg.rho_2
	pdg.UPDATE_STATE_VECTOR(ns.predPos, ns.predVelo, ns.predMass, ns.nextEta)

	tSolve         = pdg.MIN_DISTANCE()
	print("Found tSolve")
	pdg.dt, ns.dt  = .2, .2
	count, initPdg = 0, []

	while ns.startPDG:
		# calls pdg to optimize fuel use and recalculate path, tells process time
		# that there is a new solution

		t_pdgStart = ns.met

		pdg.UPDATE_STATE_VECTOR(ns.predPos, ns.predVelo, ns.predMass, ns.nextEta)

		# if tSolve is less than 6.0, replace min distance constraint with 2x current
		# distance, allows the vessel to do a vertical descent landing
		if tSolve < 6.0:
			if np.linalg.norm(pdg.x[1:3]) > pdg.dMax:
				pdg.dMax = np.linalg.norm(pdg.x[1:3])

		ns.sol, ns.eta = pdg.MIN_FUEL(tSolve=tSolve)
		print("----", tSolve, pdg.dMax)

		# tells process time that there is a new solution, decrement tSolve based on how much 
		# time has passed, and adjust predTime based on how long the prev pdg took to solve
		tSolve     -= ns.met - t_pdgStart
		ns.predTime = ns.met - t_pdgStart
		ns.new_eta  = True

		if count == 0:
			initPdg = ns.sol

		count += 1

	initPdgData = open("data/initial_pdg_data.txt", "w")
	pdgStr = ''
	for s in initPdg:
		pdgStr += str(s) + ","
	initPdgData.write(pdgStr)
	initPdgData.close()

if __name__ == '__main__':
	lock = Lock()
	manager = Manager()
	ns = manager.Namespace()

	#initiating and starting two processes
	Process_time = Process(target=process_time, name='TIME')
	Process_guid = Process(target=process_guid, name='GUID')
	Process_vess = Process(target=process_vess, name='VESS')

	Process_time.start()
	Process_guid.start()
	Process_vess.start()
	Process_time.join()
	Process_guid.join()
	Process_vess.join()