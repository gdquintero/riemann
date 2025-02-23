import numpy as np
import pandas as pd

folder = "test_alpha_0.5_1"

data1 = np.loadtxt(folder+"/output_tol_0.05.txt",delimiter=" ")
data2 = np.loadtxt(folder+"/output_tol_0.005.txt",delimiter=" ")
data3 = np.loadtxt(folder+"/output_tol_0.0005.txt",delimiter=" ")

# # Redondear todos los valores a 3 decimales
# data1 = np.round(data1, 2)
# data2 = np.round(data2, 2)
# data3 = np.round(data3, 2)

# # Convertir la tercera columna a enteros
# data1[:, 2] = data1[:, 2].astype(int)
# data2[:, 2] = data2[:, 2].astype(int)
# data3[:, 2] = data3[:, 2].astype(int)

# Aplicar np.clip en la última columna para truncar valores mayores a 3600
data1[:, -1] = np.clip(data1[:, -1], None, 3600)
data2[:, -1] = np.clip(data2[:, -1], None, 3600)
data3[:, -1] = np.clip(data3[:, -1], None, 3600)

table = np.zeros((50,18))

# Llenar la tabla con los datos en dos bloques por conjunto
table[:,0:3] = data1[:50,[0,2,3]]  # Primeras 50 filas de data1 en columnas 1-3
table[:,3:6] = data1[50:,[0,2,3]]  # Ultimas 50 filas de data1 en columnas 3-6

table[:,6:9] = data2[:50,[0,2,3]]  # Primeras 50 filas de data1 en columnas 1-3
table[:,9:12] = data3[50:,[0,2,3]]  # Ultimas 50 filas de data1 en columnas 3-6

table[:,12:15] = data3[:50,[0,2,3]]  # Primeras 50 filas de data1 en columnas 1-3
table[:,15:18] = data3[50:,[0,2,3]]  # Ultimas 50 filas de data1 en columnas 3-6

# Definir las columnas que deben guardarse como enteros
int_columns = [1, 4, 7, 10, 13, 16]

# Convertir `table` en un DataFrame de Pandas
df = pd.DataFrame(table)

# Aplicar formato correcto: enteros en columnas específicas
for col in int_columns:
    df[col] = df[col].astype(int)

# Guardar con `&` como separador (sin espacios)
df.to_csv("table.txt", sep="&", index=False, header=False, float_format="%.2f")

# Reemplazar "&" por " & " para agregar espacios
with open("table.txt", "r") as file:
    content = file.read().replace("&", " & ")

# Sobrescribir el archivo con el formato corregido
with open("table.txt", "w") as file:
    file.write(content)