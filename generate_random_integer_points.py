import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

np.random.seed()
points_20x20 = np.random.randint(0, 21, size=(20, 2))

plt.figure(figsize=(8, 8))
plt.scatter(points_20x20[:, 0], points_20x20[:, 1], c='blue', marker='o', s=100)
plt.title("20 Random Points in a 20x20 Grid")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.grid(True)
plt.xlim(-1, 21)
plt.ylim(-1, 21)

manhattan_distances_20x20 = np.zeros((20, 20))
for i in range(20):
    for j in range(20):
        manhattan_distances_20x20[i, j] = np.sum(np.abs(points_20x20[i] - points_20x20[j]))

manhattan_df_20x20 = pd.DataFrame(manhattan_distances_20x20)
manhattan_df_20x20.to_csv(r'manhattan_distances_20x20.csv', index=False)


plt.show()
