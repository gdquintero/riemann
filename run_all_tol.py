#!/usr/bin/env python3

import subprocess
import time

# Ejecutar los programas al mismo tiempo
p1 = subprocess.Popen(["python", "riemann_tol_0.5.py"])
p2 = subprocess.Popen(["python", "riemann_tol_0.05.py"])
p3 = subprocess.Popen(["python", "riemann_tol_0.005.py"])
p3 = subprocess.Popen(["python", "riemann_tol_0.0005.py"])

# Lista de procesos
processes = [p1, p2, p3]

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