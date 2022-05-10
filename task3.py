import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    pass

#os.system("/Users/Misha/Desktop/NM/Task3/a.out")

output_source = "/Users/Misha/Desktop/output"
file = open(output_source, 'r')
lines = file.readlines()
file.close()

lines = [line for line in [line.strip().split(',') for line in lines]]

# euler1
x_1  = [float(num) for num in lines[0]] # x
y_1  = [float(num) for num in lines[1]] # y(x)
y1_1 = [float(num) for num in lines[2]] # y'(x)

#euler2
x_2 = [float(num) for num in lines[3]] # x
y_2 = [float(num) for num in lines[4]] # y(x)
y1_2 = [float(num) for num in lines[5]] # y'(x)

#runge2
x_3 = [float(num) for num in lines[6]] # x
y_3 = [float(num) for num in lines[7]] # y(x)
y1_3 = [float(num) for num in lines[8]] # y'(x)

#runge4
x_4 = [float(num) for num in lines[9]] # x
y_4 = [float(num) for num in lines[10]] # y(x)
y1_4 = [float(num) for num in lines[11]] # y'(x)

#adams
x_5 = [float(num) for num in lines[12]] # x
y_5 = [float(num) for num in lines[13]] # y(x)
y1_5 = [float(num) for num in lines[14]] # y'(x)

fig1, ax1 = plt.subplots()
ax1.plot(x_1, y_1, c='green', label='euler')
ax1.plot(x_2, y_2, c='blue', label='euler_')
ax1.plot(x_3, y_3, c='red', label='runge2')
ax1.plot(x_4, y_4, c='orange', label='runge4')
ax1.plot(x_5, y_5, c='black', label='adams')
ax1.legend(loc='upper right')
ax1.grid(True)
ax1.margins(0)
plt.xlabel("x")
plt.ylabel("y(x)")
plt.gcf().set_dpi(300)
plt.savefig("y(x).png")

fig2, ax2 = plt.subplots()
ax2.plot(x_1, y1_1, c='green', label='eul')
ax2.plot(x_2, y1_2, c='blue', label='eul_recalc')
ax2.plot(x_3, y1_3, c='red', label='runge2')
ax2.plot(x_4, y1_4, c='orange', label='runge4')
ax2.plot(x_5, y1_5, c='black', label='adams')
ax2.legend(loc='upper right')
ax2.grid(True)
ax2.margins(0)
plt.xlabel("x")
plt.ylabel("y'(x)")
plt.gcf().set_dpi(300)
plt.savefig("y'(x).png")
