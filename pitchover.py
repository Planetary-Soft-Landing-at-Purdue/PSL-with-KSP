import krpc, time



def highlight(parts, color):
    if color == 'red':
        color = 1.0, 0.0, 0.0
    elif color == 'green':
        color = 0.0, 1.0, 0.0
    elif color == 'blue':
        color = 0.0, 0.0, 1.0
    for part in parts:
        part.highlight_color = color
        part.highlighted = True
    time.sleep(1.0)
    for part in parts:
        part.highlight_color = 0.0, 1.0, 0.0
        part.highlighted = False


# makes all necessary connections with ksp
conn = krpc.connect(address='10.192.30.186')
sc = conn.space_center
vessel = sc.active_vessel
ap = vessel.auto_pilot
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
g = body.surface_gravity


s1_tank = vessel.parts.with_tag('S1_Tank')[0]
rp1 = conn.add_stream(getattr, s1_tank.resources.with_resource('Kerosene')[0], 'amount')

center = vessel.parts.with_tag('0')[0]
s1_1 = vessel.parts.with_tag('1')
s1_2 = vessel.parts.with_tag('2')
s1_3 = vessel.parts.with_tag('3')
s1_4 = vessel.parts.with_tag('4')
s2_5 = vessel.parts.with_tag('5')

while rp1() > 25700: pass

center.engine.active = False
s1_1[0].engine.active = False
s1_1[1].engine.active = False
s1_2[0].engine.active = False
s1_2[1].engine.active = False
s1_3[0].engine.active = False
s1_3[1].engine.active = False
s1_4[0].engine.active = False
s1_4[1].engine.active = False

vessel.control.activate_next_stage()

center.engine.active = True

for vessel in sc.vessels:
    if len(vessel.parts.with_tag('0')) == 1: s1 = vessel
    if len(vessel.parts.with_tag('5')) == 1: s2 = vessel
s1.control.throttle = 0.5
s2.control.throttle = 0.36

y0 = s2.position(s1.reference_frame)[1]
while s2.position(s1.reference_frame)[1] < y0 + 15: pass

s2.control.throttle = 0.50
s1_1[0].engine.active = True
s1_1[1].engine.active = True

s2.auto_pilot.reference_frame = bcbf
s2.auto_pilot.target_direction = s2.direction(bcbf)
s2.auto_pilot.engage()

s1.auto_pilot.reference_frame = bcbf
s1.auto_pilot.target_direction = (
    -s1.direction(bcbf)[0],
    -s1.direction(bcbf)[1],
    -s1.direction(bcbf)[2])
s1.auto_pilot.engage()

time.sleep(3)
s1.control.throttle = 0.36

while True: pass
