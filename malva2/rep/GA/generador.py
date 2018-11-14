import sys
import random
import math
from collections import defaultdict

try:
	#Constantes
	#mayor_duracion_horas_tarea=200
	#menor_duracion_horas_tarea=12
	#mayor_disponibilidad_diaria_horas_empleado=10
	#menor_disponibilidad_diaria_horas_empleado=1
	#sueldo_hora_base_empleado_sin_experiencia=110
	#sueldo_valor_maximo_bono_por_buen_empleado_hora=23
	#deadline_dias=0
	
	#Sanity check
	if (len(sys.argv)<4):
		print "Ejecutar con la siguiente linea: ./generador <cant_vehiculos> <largo_mapa> <ancho_mapa> <min_autonomia> <max_autonomia> <cant_obstaculos> <cant_recargas>  <archivo_salida>"
		sys.exit(1)

	cant_vehiculos = int(sys.argv[1])
	largo_mapa = int(sys.argv[2])
	ancho_mapa = int(sys.argv[3])
	min_autonomia = int(sys.argv[4])
	max_autonomia = int(sys.argv[5])
	cant_obstaculos = int(sys.argv[6])
	cant_recargas = int(sys.argv[7])
	archivo_salida = sys.argv[8]

	if (cant_vehiculos<=1):
		print "La cantidad de vehiculos debe ser mayor 1."
		sys.exit(1)

	
	if (largo_mapa<=1):
		print "El largo del mapa debe ser mayor a 1."
		sys.exit(1)

	if (ancho_mapa<=1):
		print "El ancho del mapa debe ser mayor a 1."
		sys.exit(1)


	#VARIABLES:
	
	#Autonomia de los vehiculos, toma valores comprendidos entre min y max pasados por parametro
	autonomia_vehiculo = []	

	#Se define el mapa tomando los siguintes valores:
	# 1 visitable
	# 0 no visitable
	# 2 punto de carga
	mapa = [[0 for x in range(largo_mapa)] for y in range(ancho_mapa)]	
	
	for i in range(0,cant_vehiculos):
		autonomia_vehiculo.append(random.randint(min_autonomia,max_autonomia))


	for i in range(0,largo_mapa):
		for j in range(0,ancho_mapa):
			mapa[j][i] = 1
	
	contrec = 0
	contobs = 0				
	while (contobs < cant_obstaculos) or (contrec < cant_recargas):
		i = random.randint(0,largo_mapa-1)
		j = random.randint(0,ancho_mapa-1)
		aux = random.randint(0,2)
		if (aux == 2):
			if (contrec < cant_recargas):
				mapa[j][i] = 2
				contrec = contrec + 1
			else:
				aux = random.randint(0,1)
				if (aux == 0):
					if (contobs > cant_obstaculos):
						mapa[j][i] = 0 
						contobs = contobs + 1
					else:
						mapa[j][i] = 1
				else:
					mapa[j][i] = 1		
		elif (aux == 0):
			if (contobs < cant_obstaculos):
				mapa[j][i] = 0 
				contobs = contobs + 1
			else:
				aux = random.randint(1,2)
				if (aux == 2):
					if (contrec < cant_recargas):
						mapa[j][i] = 2
						contrec = contrec + 1	
					else:	
						mapa[j][i] = 1	
				else:
					mapa[j][i] = 1								

	
	#archivo de salida de empleados
	archivo = open ("%s"%archivo_salida, "w")
	archivo.write("%d"%cant_vehiculos)
	archivo.write("\n")
	archivo.write("%d"%largo_mapa)
	archivo.write("\n")
	archivo.write("%d"%ancho_mapa)
	archivo.write("\n")
	
	archivo.write("%d"%autonomia_vehiculo[0])
	for i in range(1,cant_vehiculos):
		archivo.write(" %d"%autonomia_vehiculo[i])
	archivo.write("\n")
	
	for i in range(0,largo_mapa):
		archivo.write("%d"%mapa[0][i])
		for j in range(0,ancho_mapa):
			archivo.write(" %d"%mapa[j][i])	
		archivo.write("\n")		
	
	

except IOError as error:
	print error
