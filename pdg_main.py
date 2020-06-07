from multiprocessing import Process,Value,Lock
from config import *
import krpc,time


def process_time():
	#  connect to server, set conn & vessel
	#  create & start MET stream, run update_met on callback
	#  wait for termination
	try:
		conn = krpc.connect(address='192.168.1.169')
		vessel = conn.space_center.active_vessel
		met_stream = conn.add_stream(getattr,vessel,'met')
		met_stream.add_callback(update_met)
		met_stream.start()
		while True: pass
	#  if connection fails, simulate MET updates
	except:
		print('Connection failed, simulating MET')
		while True:
			update_met(met.value)
			time.sleep(0.04)
			met.value += 0.04

def update_met(self):
	#  update MET and runs callback-synced functions
	met.value = self
	control()

def control():
	#  send control commands
	#  if new solution available, reset control elapsed time
	global cet0
	if new_eta.value == 1:
		cet0 = met.value
		cet = 0
		new_eta.value = 0
	cet = met.value-cet0
	with lock:
		print('--- cet:',round(cet,3))


def process_guidance():
	#  initialize guidance, call pdg loop while T > 0
	#  terminate process_time when finished
	T = 1
	while met.value == 0: pass
	while T > 0:
		t0 = met.value
		time.sleep(0.25)
		dT = met.value-t0
		with lock:
			print('GdT:',round(dT,3),';','T:',round(T,2))
			eta = [T,T]
			print(eta[:])
		T -= dT
		new_eta.value = 1
	process_time.terminate()


if __name__ == '__main__':

	lock = Lock()

	process_time = Process(target=process_time,args=(),name='Time')
	process_guidance = Process(target=process_guidance,args=(),name='Guidance')

	process_time.start()
	process_guidance.start()