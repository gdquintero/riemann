import subprocess

# Ejecutar dos programas al mismo tiempo
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

# Esperar a que ambos terminen (opcional)
p1.wait()
p2.wait()
p3.wait()
p4.wait()
p5.wait()
p6.wait()
p7.wait()
p8.wait()
p9.wait()
p10.wait()

print("Todos las instancias han terminado.")