#!/usr/bin/env python3

import subprocess
import time

# Ejecutar los programas al mismo tiempo
# p1 = subprocess.Popen(["python", "riemann_0.5.py"])
p2 = subprocess.Popen(["python", "riemann_0.05.py"])
p3 = subprocess.Popen(["python", "riemann_0.005.py"])
p4 = subprocess.Popen(["python", "riemann_0.0005.py"])