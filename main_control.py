from multiprocessing import Process, Manager, Lock
from math import log
import krpc, time, pdg, math
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

	while position_stream()[1] < 600: pass

	vessel.control.throttle = 0
	vessel.auto_pilot.sas 	= False
	vessel.auto_pilot.engage()
	vessel.auto_pilot.reference_frame = pcpf

	print("Starting controlled descent")
	startTime, n, tPrev, ns.eta = met_stream(), 0, -ns.dt, [0, 0, 0, 0]

	# continuously updates the vessel's thrust until it reaches a certain
	# altitude above its target landing spot
	while ns.new_eta == None or position_stream()[1] > 1:
		# updates ns namespace variables
		ns.position  = position_stream()
		ns.stateVect = [position_stream()[1], position_stream()[0], position_stream()[2],
						velocity_stream()[1], velocity_stream()[0], velocity_stream()[2],
						log(mass_stream()),
						ns.eta[n+0], ns.eta[n+1], ns.eta[n+2], ns.eta[n+3]
					   ]
		ns.met      = met_stream()

		ns.startPDG = True

		# waits until a new eta has been found
		if ns.new_eta == None: pass
		else:
			# if Guidance found a new path, reset startTime, n, and tPrev
			if ns.new_eta == True:
				ns.new_eta = False
				startTime, n, tPrev = met_stream(), 0, -ns.dt

			# if an entire time step has passed, move to next set of eta values
			if int((time.time() - startTime) / ns.dt) != tPrev: 
				n, tPrev = n + 4, int((met_stream() - startTime) / ns.dt)

			vessel.auto_pilot.target_direction = ns.eta[n+1], ns.eta[n], ns.eta[n+2]
			vessel.control.throttle 		   = ns.eta[n+3] * (mass_stream() / pdg.rho_2)

	print("Finished Guidance")
	ns.startPDG = False
	vessel.control.throttle = 0

def process_guid():
	global tWait, tSolve, dMax, dt

	ns.new_eta = None
	ns.dt      = .2
	# waits until process time tells it to start looking for a solution
	while ns.startPDG == False: pass

	# timeWait = pdg.PDG(ns.dt, ns.stateVect, initialSearch=True)
	# print("Will wait", timeWait, "seconds before starting pdg")
	# time.sleep(timeWait)
	# print("Stopped waiting")

	ns.dt = .2
	count = 0
	tSolve, dMax, _ = pdg.PDG(ns.dt, ns.stateVect, minDistance=True)

	while ns.startPDG:
		# for every 10 pdg calls, re-optimize descent time, don't do this if the 
		# vessel's altitude is less than 50 meters, optimize for final distance
		if count % 10 == 0 and ns.position[1] > 50:
			tSolve, dMax, _ = pdg.PDG(ns.dt, ns.stateVect, minDistance=True)
			print("----------------------------------")
			print("----------Found new path----------")
			print("----------------------------------")

		# calls pdg to optimize fuel use and recalculate path, tells process time
		# that there is a new solution
		startTime  = ns.met
		ns.eta, _  = pdg.PDG(ns.dt, ns.stateVect, tSolve=tSolve, tSolveTotal=tSolve, dMax=dMax)
		print("----", tSolve)
		ns.new_eta = True

		count  += 1
		tSolve -= ns.met - startTime
	
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


