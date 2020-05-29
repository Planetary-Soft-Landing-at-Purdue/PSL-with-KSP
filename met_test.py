import krpc
import time

#  connect to server
conn = krpc.connect(address = '192.168.1.169')

#  set active vessel
vessel = conn.space_center.active_vessel

#  create UT stream
ut = conn.add_stream(getattr, conn.space_center, 'ut')

def update_met(ut):

	#  initialize MET
	if 'met' not in globals():
		global met
		global ut_last
		met = 0
		ut_last = ut
	#  update MET
	else:
		met = met + ut - ut_last
		ut_last = ut

#  create and start MET callback
ut.add_callback(update_met)
ut.start()

