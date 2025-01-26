import subprocess

# Crear una lista de nombres de programas
programas = [f"riemann_{i}.py" for i in range(1, 21)]  # Genera riemann_1.py, ..., riemann_20.py

# Crear procesos para cada programa
procesos = [subprocess.Popen(["python", programa]) for programa in programas]

# Esperar a que todos los procesos terminen
for p in procesos:
    p.wait()

print("Todas las instancias han terminado.")