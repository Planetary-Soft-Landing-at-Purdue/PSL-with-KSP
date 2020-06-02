from multiprocessing import Process, Value, Lock
import time, krpc


#  process aligned with timesteps
#	connect to server
#	set conn and vessel
#	create and start MET stream 
#	on callback, run update_met()
#	wait for termination (from Process_guidance)
def Process_time():
	global met_stream
	conn = krpc.connect(address='192.168.1.169')
	vessel = conn.space_center.active_vessel
	met_stream = conn.add_stream(getattr,vessel,'met')
	met_stream.add_callback(update_met)
	met_stream.start()
	while True:
	 	pass

#  update shared Value met, global to Process_time
def update_met(self):
	global met_stream
	met.value = met_stream()

	#  everything to be run on each timestep
	control()

#  guidance loop
def Process_guidance():

	#  initialize guidance variables
	global new_eta
	T = 1

	while met.value == 0: pass

	#  call pdg loop
	while T > 0:
		t0 = met.value
		time.sleep(0.25)
		dT = met.value-t0
		lock.acquire()
		print('GdT:',round(dT,3),';','T:',round(T,2))
		lock.release()
		T -= dT
		new_eta.value = 1
	Process_time.terminate()

#  control function, run on each timestep
#	if new solution available from Process_guidance(),
#	then reset control elapsed time (cet)
#	print cet: placeholder for sending control command
def control():
	global new_eta, cet0
	if new_eta.value == 1:
		cet0 = met.value
		cet = 0
		new_eta.value = 0
	cet = met.value-cet0
	lock.acquire()
	print('--- cet:',round(cet,3))
	lock.release()


if __name__ == "__main__":

	lock = Lock()
	met = Value('f')
	new_eta = Value('i')
	new_eta.value = 1

	Process_time = Process(target=Process_time, args=(), name='Time')
	Process_guidance = Process(target=Process_guidance, args=(), name='Guidance')

	Process_time.start()
	Process_guidance.start()