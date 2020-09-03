import krpc, time

def is_opposite_engine(engine1, engine2):
	vessel = engine1.part.vessel
	rf = vessel.reference_frame
	pos1 = engine1.part.position(rf)
	pos2 = engine2.part.position(rf)

	if abs(pos1[0]+pos2[0]) < 0.01 and abs(pos1[2]+pos2[2]) < 0.01:
		highlight((engine1.part, engine2.part), 'green')
		return True
	else:
		highlight((engine1.part, engine2.part), 'red')
		return False

def highlight(parts, color):
	if color == 'red': color = 1.0, 0.0, 0.0
	elif color == 'green': color = 0.0, 1.0, 0.0
	elif color == 'blue': color = 0.0, 0.0, 1.0
	for part in parts:
		part.highlight_color = color
		part.highlighted = True
	time.sleep(0.05)
	for part in parts:
		part.highlight_color = 0.0, 1.0, 0.0
		part.highlighted = False

# makes all necessary connections with ksp
conn = krpc.connect(address='10.192.28.146')
sc = conn.space_center
vessel = sc.active_vessel
ap = vessel.auto_pilot
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
g = body.surface_gravity

for vessel in sc.vessels:
  for engine in vessel.parts.engines:
      if engine.part.name == "ROE-Merlin1D":
          s1 = vessel
      elif engine.part.name == "ROE-Merlin1DV":
          s2 = vessel

s2.control.throttle = 1
s2.auto_pilot.engage()
s2.auto_pilot.reference_frame = bcbf
s2.auto_pilot.target_direction = s2.direction(bcbf)

s1.control.throttle = 0
center = 0; left = 0; right = 0
s1rf = s1.reference_frame
for engine in s1.parts.engines:
	pos = engine.part.position(s1rf)
	if center == 0:
		if abs(pos[0]) < 0.001 and abs(pos[2]) < 0.001:
			center = engine
			print(center.part.position(s1rf))
	elif left == 0: left = engine
	elif is_opposite_engine(left, engine): right = engine