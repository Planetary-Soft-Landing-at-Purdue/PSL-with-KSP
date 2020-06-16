import krpc
from numpy import arctan2
from math import sqrt

#  BCBF position tuple to coordinate tuple
def pos_to_coords(pos):
	lat = arctan2(pos[1],pos[0])
	lon = arctan2(pos[2],pos[0])
	return lat,lon

#  find impact time of suborbital vessel
def find_t_impact_g(vessel,space_center):
	orbit = vessel.orbit
	body = orbit.body
	bcbf = body.reference_frame
	t0 = space_center.ut
	t_impact = space_center.ut
	x_impact = vessel.orbit.position_at(t_impact,bcbf)
	c_impact = pos_to_coords(x_impact)

	#  increment impact time until impact point found
	while body.altitude_at_position(orbit.position_at(t_impact,bcbf),bcbf) > 0:
		t_impact += 10
		x_impact = vessel.orbit.position_at(t_impact,bcbf)
		c_impact = pos_to_coords(x_impact)
		if t_impact-t0 > orbit.period:
			print('Error in find_t_impact: No impact found')
			return
	return t_impact-t0 

def zero_quadratic(a,b,c):
	try:
		return (
			(-b+sqrt(b**2-4*a*c))/(2*a),
			(-b-sqrt(b**2-4*a*c))/(2*a))
	except:
		return 'No roots found'

#  Schedules optimal ignition time for PDG
def schedule_pdg():
	t_ignition = t_current
	while True:
		state0 = state(t_ignition)					#  state at this ignition time
		solvetime = goldenSearch(state0,t_ignition)	#  optimal solvetime for this ignition time
		fuel = pdg(state0,solvetime)[0]				#  fuel use for optimal solvetime with this ignition time
		if fuel_last not in globals or fuel < fuel_last:
			t_ignition += 1
		else:
			return t_ignition

#  test code here
conn = krpc.connect(address='192.168.1.169')
vessel = conn.space_center.active_vessel
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
x0 = vessel.flight(bcbf).mean_altitude
v0 = vessel.flight(bcbf).vertical_speed
g = body.surface_gravity
print(zero_quadratic(-g/2,v0,x0)[1])