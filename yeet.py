from multiprocessing import Value
import krpc, time, datetime


def met_callback(self):
	met.value = self


def position_callback(self):
	x.value = self[0]
	y.value = self[1]
	z.value = self[2]


conn = krpc.connect(address='192.168.1.181')
sc = conn.space_center
vessel = sc.active_vessel
orbit = vessel.orbit
body = orbit.body
bcbf = body.reference_frame
pcpr = conn.space_center.ReferenceFrame.create_relative(
	bcbf,
	position=vessel.position(bcbf),
	rotation=vessel.rotation(bcbf))

x = Value('d')
y = Value('d')
z = Value('d')
met = Value('d')
x.value = 0
y.value = 0
z.value = 0
met.value = 0

position = conn.add_stream(vessel.position,pcpr)
position.add_callback(position_callback)
position.start()

met = conn.add_stream(getattr, vessel, 'met')
met.add_callback(met_callback)
met.start()

while True:
	print(
		'{0: 3f}'.format(x.value),
		'{0: 3f}'.format(y.value),
		'{0: 3f}'.format(z.value),
		' ',
		str(datetime.timedelta(seconds=met.value))
		)
	time.sleep(1)