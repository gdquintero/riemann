#!/usr/bin/env python3

import subprocess
import time

# Ejecutar los programas al mismo tiempo
p1 = subprocess.Popen(["python", "riemann_0.01.py"])
p2 = subprocess.Popen(["python", "riemann_0.001.py"])
p3 = subprocess.Popen(["python", "riemann_0.0001.py"])

# Esperar a que todos terminen
p1.wait()
p2.wait()
p3.wait()

print("Todos los programas terminaron.")