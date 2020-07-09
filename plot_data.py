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

control = open("pdg_data.txt", "r")
control = control.readlines()[0].split("[")[1]
control = control.split("]")[0]
controlData = control.split(",")

real = open("real_data.txt", "r")
real = real.readlines()[0].split("[")[1]
real = real.split("]")[0]
realData = real.split(",")

init = open("initial_pdg_data.txt", "r")
init = init.readlines()[0]
initData = init.split(",")[11:]

t_listControl = range(len(controlData) // 11 - 1)
t_listReal = range(len(realData) // 11 - 1)
t_listInit = range(len(initData) // 11)

c_z, c_x, c_y, c_dz, c_dx, c_dy, c_z_t, c_tz, c_tx, c_ty, c_s = seperateData(controlData)
r_z, r_x, r_y, r_dz, r_dx, r_dy, r_z_t, r_tz, r_tx, r_ty, r_s = seperateData(realData)
i_z, i_x, i_y, i_dz, i_dx, i_dy, i_z_t, i_tz, i_tx, i_ty, i_s = seperateData(initData)

plt.figure(1)
plt.plot(t_listControl, c_z, t_listReal, r_z, t_listInit, i_z)
plt.title("Vertical Displacement over time")
plt.savefig('figures/z_displacement')

plt.figure(2)
plt.plot(t_listControl, c_tz, t_listReal, r_tz, t_listInit, i_tz)
plt.title("z_thrust over time")
plt.savefig('figures/z_thrust')

plt.figure(3)
plt.plot(t_listControl, c_dz, t_listReal, r_dz, t_listInit, i_dz)
plt.title("velocity over time")
plt.savefig('figures/z_velocity')

plt.figure(4)
plt.plot(t_listControl, c_z_t, t_listReal, r_z_t, t_listInit, i_z_t)
plt.title("mass over time")
plt.savefig('figures/mass')

plt.figure(5)
plt.plot(t_listControl, c_s, t_listReal, r_s, t_listInit, i_s)
plt.title("sigma over time")
plt.savefig('figures/sigma')

plt.figure(6)
plt.plot(t_listControl, c_x, t_listReal, r_x,  t_listInit, i_x)
plt.title("X Displacement over time")
plt.savefig('figures/x_displacement')

plt.figure(7)
plt.plot(t_listControl, c_y, t_listReal, r_y, t_listInit, i_y)
plt.title("Y Displacement over time")
plt.savefig('figures/y_displacement')

plt.figure(8)
plt.plot(t_listControl, c_tx, t_listReal, r_tx,  t_listInit, i_tx)
plt.title("x_thrust over time")
plt.savefig('figures/x_thrust')

plt.figure(9)
plt.plot(t_listControl, c_ty, t_listReal, r_ty, t_listInit, i_ty)
plt.title("y_thrust over time")
plt.savefig('figures/y_thrust')

plt.figure(10)
plt.plot(t_listControl, c_dx, t_listReal, r_dx, t_listInit, i_dx)
plt.title("x velocity over time")
plt.savefig('figures/x_velocity')

plt.figure(11)
plt.plot(t_listControl, c_dy, t_listReal, r_dy, t_listInit, i_dy)
plt.title("y velocity over time")
plt.savefig('figures/y_velocity')
