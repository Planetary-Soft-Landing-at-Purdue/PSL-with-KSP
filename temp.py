import krpc, time


def is_opposite_engine(engine1, engine2):
    vessel = engine1.part.vessel
    rf = vessel.reference_frame
    pos1 = engine1.part.position(rf)
    pos2 = engine2.part.position(rf)

    if abs(pos1[0] + pos2[0]) < 0.01 and abs(pos1[2] + pos2[2]) < 0.01:
        return True
    else:
        return False


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
conn = krpc.connect()
sc = conn.space_center
vessel = sc.active_vessel
ap = vessel.auto_pilot
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
g = body.surface_gravity

rp1 = conn.add_stream(vessel.resources_in_decouple_stage(0).amount, 'Kerosene')

while rp1() > 27567: pass

vessel.control.throttle = 0.0
vessel.control.activate_next_stage()

for vessel in sc.vessels:
    for engine in vessel.parts.engines:
        if engine.part.name == "ROE-Merlin1D":
            s1 = vessel
        elif engine.part.name == "ROE-Merlin1DV":
            s2 = vessel

s2.control.throttle = 1.0
s2.auto_pilot.engage()
s2.auto_pilot.reference_frame = bcbf
print(s2.direction(bcbf))
s2.auto_pilot.target_direction = s2.direction(bcbf)

s1.control.throttle = 0.0
center = 0; left = 0; right = 0
s1rf = s1.reference_frame
for engine in s1.parts.engines:
    pos = engine.part.position(s1rf)
    if center == 0:
        if abs(pos[0]) < 0.001 and abs(pos[2]) < 0.001:
            center = engine
    elif left == 0:
        left = engine
    elif is_opposite_engine(left, engine):
        right = engine
highlight((center.part,), 'Red')
highlight((left.part, right.part), 'Blue')
