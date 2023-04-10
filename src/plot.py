import matplotlib.pyplot as plt

def generate_graph(data, timepoints):
    plt.figure(figsize=(14, 9))
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

def generate_extinction_graph(simulated_data, true_data):
    site_count = len(simulated_data) - 1
    plt.figure(figsize=(15, 10))
    plt.grid(visible = True)
    plt.xlim(0, site_count)
    plt.ylim(0, 1)
    plt.xlabel('number of methylated sites that cell has',fontsize=20)
    plt.ylabel('extinction probability starting with a single cell',fontsize=20)

    plt.plot(range(site_count + 1), simulated_data, label = 'monte carlo simulated extinction rate')
    plt.plot(range(site_count + 1), true_data, label = 'calculated extinction rate')
    plt.legend(loc = 'lower left')

    plt.show()

