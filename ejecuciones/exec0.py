import os
import subprocess


script_ejecuciones = "/home/giovani/AE/ejecuciones/script_ejecuciones"
prefixInstance = "instancias/instancia_"
prefixNewGA = "ArchivosCFG/"
posfixNewGA = "/newGA.cfg"
for j in range(1, 28):
    newGA = prefixNewGA + str(j) + posfixNewGA
    for i in range(1, 6):
        instance = prefixInstance + str(i)
        for k in range(4):
            args = (script_ejecuciones, instance, newGA)
            # print ("instance: ",instance," newGA: ", newGA)
            popen = subprocess.Popen(args, stdout=subprocess.PIPE)
            popen.wait()
            output = popen.stdout.read()
            # print(output)
            # print("terminaEjecucion")
            # os.system('script_ejecuciones + " " + instance + " " + newGA)
