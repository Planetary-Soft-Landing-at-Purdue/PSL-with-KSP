from multiprocessing import Value

met = Value('d')
met.value = 0

mass = Value('d')
mass.value = 0

new_eta = Value('i')
new_eta.value = 1

position = Value('d'),Value('d'),Value('d')
velocity = Value('d'),Value('d'),Value('d')