import krpc
from numpy import arctan2
from numpy.linalg import norm

#  BCBF position tuple to coordinate tuple
def pos_to_coords(pos):
	lat = arctan2(pos[1],pos[0])
	lon = arctan2(pos[2],pos[0])
	return lat,lon

#  find impact time of suborbital vessel
def find_t_impact(vessel,space_center):
	orbit = vessel.orbit
	body = orbit.body
	bcbf = body.reference_frame
	t0 = space_center.ut
	t_impact = space_center.ut
	x_impact = vessel.orbit.position_at(t_impact,bcbf)
	c_impact = pos_to_coords(x_impact)

	#  increment impact time until impact point found
	while (orbit.radius_at(t_impact) >
		norm(body.surface_position(c_impact[0],c_impact[1],bcbf))):
		t_impact += 10
		x_impact = vessel.orbit.position_at(t_impact,bcbf)
		c_impact = pos_to_coords(x_impact)
		if t_impact-t0 > orbit.period:
			print('Error in find_t_impact: No impact found')
			return
	return t_impact-t0


#  test code here
conn = krpc.connect(address='192.168.1.169')
vessel = conn.space_center.active_vessel
print('Time to impact: ',find_t_impact(vessel,conn.space_center))