from multiprocessing import Process, Value
import time

#  MET stream simulator
def Process_time():

	#  MET change callback
	#  do everything for this timestep
	while True:
		time.sleep(0.1)
		print('MET:',round(met.value, 3),'; CET:',control())
		met.value += 0.1

#  guidance loop
def Process_guidance():

	#  initialize guidance variables
	global new_eta
	T = 3

	#  call pdg loop
	while T > 0:
		t0 = met.value
		time.sleep(0.5)
		dT = met.value-t0
		print('GdT:',round(dT,3),';','T:',round(T,3))
		T -= dT
		new_eta.value = 1
	Process_time.terminate()

def control():
	
	global new_eta, cet0
	if new_eta.value == 1:
		cet0 = met.value
		cet = 0
		new_eta.value = 0
	cet = met.value-cet0
	return(round(cet,3))


if __name__ == "__main__":

	met = Value('f')
	new_eta = Value('i')
	new_eta.value = 1

	Process_time = Process(target=Process_time, args=(), name='Time')
	Process_guidance = Process(target=Process_guidance, args=(), name='Guidance')

	Process_time.start()
	Process_guidance.start()