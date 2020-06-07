from multiprocessing import Value,Array

met = Value('d')
met.value = 0

new_eta = Value('i')
new_eta.value = 1

eta = Array('d', 0)