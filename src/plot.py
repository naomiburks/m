import matplotlib.pyplot as plt

def generate_graph(data, timepoints):
    plt.figure(figsize=(15, 10))
    plt.grid(visible = True)
    plt.xlim(0, timepoints[-1])
    plt.ylim(0, max([max(row) for row in data]))
    plt.xlabel('time',fontsize=20)
    plt.ylabel('cell counts',fontsize=20)

    
    for i in range(len(data[0])):
        pop_counts = data[:, i]
        plt.plot(timepoints, pop_counts, label = f'{i} methylated sites')

    plt.legend(loc = 'upper right')

    plt.show()

