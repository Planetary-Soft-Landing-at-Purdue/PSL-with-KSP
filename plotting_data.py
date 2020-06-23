import matplotlib.pyplot as plt

def seperateData(data):
	z, x, y, dz, dx, dy, z_t, tz, tx, ty, s = [], [], [], [], [], [], [], [], [], [], []

	for t in range(0, len(data) - 11, 11):
		z.append( 	float(data[t]) 	   )
		x.append( 	float(data[t + 1]) )
		y.append( 	float(data[t + 2]) )
		dz.append( 	float(data[t + 3]) )
		dx.append( 	float(data[t + 4]) )
		dy.append( 	float(data[t + 5]) )
		z_t.append( float(data[t + 6]) )
		tz.append( 	float(data[t + 7]) )
		tx.append( 	float(data[t + 8]) )
		ty.append( 	float(data[t + 9]) )
		s.append( 	float(data[t + 10]))

	return z, x, y, dz, dx, dy, z_t, tz, tx, ty, s

control = open("dataFile_control.txt", "r")
controlData = control.readlines()[0].split(",")
control.close()

real = open("dataFile_real.txt", "r")
realData = real.readlines()[0].split(",")
real.close()

t_listControl = range(len(controlData) // 11)
t_listReal = range(len(realData) // 11)

c_z, c_x, c_y, c_dz, c_dx, c_dy, c_z_t, c_tz, c_tx, c_ty, c_s = seperateData(controlData)
r_z, r_x, r_y, r_dz, r_dx, r_dy, r_z_t, r_tz, r_tx, r_ty, r_s = seperateData(realData)

plt.figure(1)
plt.plot(t_listControl, c_z, t_listReal, r_z)
plt.title("Displacement over time")
plt.savefig('figures/displacement')

plt.figure(2)
plt.plot(t_listControl, c_tz, t_listReal, r_tz)
plt.title("thrust over time")
plt.savefig('figures/thrust')

plt.figure(3)
plt.plot(t_listControl, c_dz, t_listReal, r_dz)
plt.title("velocity over time")
plt.savefig('figures/velocity')

plt.figure(4)
plt.plot(t_listControl, c_z_t, t_listReal, r_z_t)
plt.title("mass over time")
plt.savefig('figures/mass')

plt.figure(5)
plt.plot(t_listControl, c_s, t_listReal, r_s)
plt.title("sigma over time")
plt.savefig('figures/sigma')
