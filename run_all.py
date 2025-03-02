#!/usr/bin/env python3

import subprocess
import time

# Ejecutar los programas al mismo tiempo
p1 = subprocess.Popen(["python", "riemann_1.py"])
p2 = subprocess.Popen(["python", "riemann_2.py"])
p3 = subprocess.Popen(["python", "riemann_3.py"])
p4 = subprocess.Popen(["python", "riemann_4.py"])
p5 = subprocess.Popen(["python", "riemann_5.py"])
p6 = subprocess.Popen(["python", "riemann_6.py"])
p7 = subprocess.Popen(["python", "riemann_7.py"])
p8 = subprocess.Popen(["python", "riemann_8.py"])
p9 = subprocess.Popen(["python", "riemann_9.py"])
p10 = subprocess.Popen(["python", "riemann_10.py"])
p11 = subprocess.Popen(["python", "riemann_11.py"])
p12 = subprocess.Popen(["python", "riemann_12.py"])
p13 = subprocess.Popen(["python", "riemann_13.py"])
p14 = subprocess.Popen(["python", "riemann_14.py"])
p15 = subprocess.Popen(["python", "riemann_15.py"])
p16 = subprocess.Popen(["python", "riemann_16.py"])
p17 = subprocess.Popen(["python", "riemann_17.py"])
p18 = subprocess.Popen(["python", "riemann_18.py"])
p19 = subprocess.Popen(["python", "riemann_19.py"])
p20 = subprocess.Popen(["python", "riemann_20.py"])

# Lista de procesos
processes = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,20]

try:
    while True:
        # Verificar si algún proceso sigue en ejecución
        running = any(process.poll() is None for process in processes)
        if not running:
            break  # Salir del bucle si todos los procesos han terminado

        # Mostrar mensaje de "Rodando"
        print("Rodando...", end="\r", flush=True)
        time.sleep(1)  # Esperar 1 segundo antes de verificar nuevamente
except KeyboardInterrupt:
    print("\nInterrumpido por el usuario. Terminando procesos...")
    for process in processes:
        process.terminate()

# Confirmar que han terminado
print("Todos los programas han finalizado.")