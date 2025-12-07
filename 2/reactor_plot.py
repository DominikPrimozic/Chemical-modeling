
import subprocess
import numpy as np
import matplotlib.pyplot as plt


file_path = "input/reactor.txt"
params = {}


with open(file_path, 'r') as file:
    params['k'], params['P'], params['s'], params['C'], params['T0'] = map(float, file.readline().split())
    params['x0']= float(file.readline())
    params['steps'] = int(file.readline())
    params['dt'] = float(file.readline())
    params['file_output'] = file.readline().strip()


output_file = f"output/reactor/{params['file_output']}.txt"


subprocess.run(["prva_naloga.exe"], check=True)



t, x = np.loadtxt(output_file , unpack=True)

plt.figure(figsize=(10, 5))
plt.plot(t, x, label="x(t)")
plt.xlabel("Čas")
plt.ylabel("Temperatura")
plt.title(f"$k$={params['k']}, $P$={params['P']}, $\\sigma$={params['s']}, $C$={params['C']}, $T_0$={params['T0']}")

plt.text(1.05, 0.5, f"$x(0)$={params['x0']}\n $\\Delta t$={params['dt']}",
         horizontalalignment='left', verticalalignment='center',
         transform=plt.gca().transAxes, fontsize=12)

plt.subplots_adjust(right=0.85)

output_image = f"output/reactor/{params['file_output']}.png"
plt.savefig(output_image, dpi=300) 
plt.show()