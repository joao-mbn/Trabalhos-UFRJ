import matplotlib.pyplot as plt

temperatures = [45,60,75,60,45,50,68,84,70,80,90,100,90,90,90,90]
cat_amount = [0.1176,0.1176,0.1176,0.1078,0.01671,0.147,0.147,0.147,0.75,0.75,0.75,0.75,0.25,1.5,0.75,0.75]
molar_ratio = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5,2]

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(temperatures, cat_amount, molar_ratio)
ax.set_xlabel('temperture (ºC)')
ax.set_ylabel('amount of catalyst (eqg H2SO4/g sol. (%)')
ax.set_zlabel('IbOH - IbAc molar ratio')
plt.show()
