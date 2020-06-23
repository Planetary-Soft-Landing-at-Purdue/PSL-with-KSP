from multiprocessing import Process, Manager, Lock
import krpc, time, pdg
import numpy as np

# initializes everything

conn 	= krpc.connect(address='192.168.0.109')
sc 		= conn.space_center
vessel 	= sc.active_vessel
orbit 	= vessel.orbit
body 	= orbit.body
bcbf 	= body.reference_frame
bci 	= body.non_rotating_reference_frame
omega 	= body.angular_velocity(bci)

pcpf = sc.ReferenceFrame.create_relative(bcbf, 
		position=vessel.position(bcbf),
		rotation=vessel.rotation(bcbf))

# starts all of the necessary connections and streams

position_stream = conn.add_stream(vessel.position, pcpf)
position_stream.start()

velocity_stream = conn.add_stream(vessel.velocity, pcpf)
velocity_stream.start()

mass_stream = conn.add_stream(getattr, vessel, 'mass')
mass_stream.start()

met_stream = conn.add_stream(getattr, vessel, 'met')
met_stream.start()

print("Started all streams")

# Activates sas and activates the launch. Contiunes
# ascending until it reaches an altitude of 400 meters, 
# then begins its descent

print("Starting Launch")

vessel.auto_pilot.sas = True
vessel.control.throttle = .66
vessel.control.activate_next_stage()

while vessel.flight().surface_altitude < 600: pass

vessel.control.throttle = 0
vessel.auto_pilot.sas 	= False

# find the current state vector, calls pdg, finds a solution

stateVect = np.array([
	 position_stream()[1], position_stream()[0], position_stream()[2],
	-velocity_stream()[1], velocity_stream()[0], velocity_stream()[2],
	np.log(mass_stream()), 0, 0, 0, 0])

dt = 1

print("Looking for solution")
tWait, tSolve, dMax = pdg.PDG(dt, stateVect, initialSearch=True)
eta, sol = pdg.PDG(dt, stateVect, tWait=tWait, tSolve=tSolve, dMax=dMax)
print("Found solution")

'''
	This is where eta is used to guide the vessel's descent. It waits
	until it reaches the altiude where it started the pdg call. Then, 
	it continues updating the vessel's direction and thrust magnitude
	until it reaches 5 meters above the ground. For every dt, it updates
	the vessel's thrust to a new vector from eta.
'''

realSol = []

vessel.auto_pilot.engage()
vessel.auto_pilot.reference_frame = pcpf

while vessel.flight().surface_altitude > 600: pass
print("Starting controlled descent")
startTime, n, tPrev = time.time(), 0, 0

while vessel.flight().surface_altitude > 5 and n < len(eta) - 4:
	# waits until a there's a new discrete time step before updating the
	# vessel's thrust
	while int( (time.time() - startTime) / dt ) == tPrev: pass

	n, tPrev, mass = n + 4, int( (time.time() - startTime) / dt ), mass_stream()

	vessel.auto_pilot.target_direction = eta[n+1], eta[n], eta[n+2]
	vessel.control.throttle 		   = eta[n+3] * (mass / pdg.rho_2)

	pos  = position_stream()
	velo = velocity_stream()
	realSol.extend([pos[1], pos[0], pos[2], velo[1], velo[0], velo[2], np.log(mass), eta[n]*eta[n+3], eta[n+1]*eta[n+3], eta[n+2]*eta[n+3], eta[n+3]])

vessel.control.throttle = 0

# here we write two files, one containing the control eta 
# (what its supposed to do), and the other contains the real
# eta (what it actually does)

pdgDataFile = open("dataFile_control.txt", 'w')
pdgString = ""
for s in sol:
	pdgString += str(s) + ','
pdgDataFile.write(pdgString)
pdgDataFile.close()

realDataFile = open("dataFile_real.txt", "w")
realString = ""
for s in realSol:
	realString += str(s) + ','
realDataFile.write(realString)
realDataFile.close()