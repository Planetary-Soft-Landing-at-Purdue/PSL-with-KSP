from multiprocessing import Value

met = Value('d')
met.value = 0

new_eta = Value('i')
new_eta.value = 1