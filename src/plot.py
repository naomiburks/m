import matplotlib.pyplot as plt

def plot_simulation(data, timepoints):
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

def plot_extinction_monte_carlo_graph(simulated_data, true_data, simulations_per_type = 100):
    site_count = len(simulated_data) - 1
    plt.figure(figsize=(15, 10))
    plt.grid(visible = True)
    plt.xlim(0, site_count)
    plt.ylim(0, 1)
    plt.xlabel('number of methylated sites that cell has',fontsize=20)
    plt.ylabel('extinction probability starting with a single cell',fontsize=20)

    true_extinction_rate_variances = [(r * (1 - r) / simulations_per_type) ** 0.5 
                                     for r in true_data]

    plt.errorbar(range(site_count + 1), true_data, 
                 yerr = true_extinction_rate_variances, 
                 label = 'calculated extinction rate',
                 capsize=6)
    plt.plot(range(site_count + 1), simulated_data, label = 'monte carlo simulated extinction rate')
    plt.legend(loc = 'lower left')

    plt.show()

def plot_extinction_limit_graph(data, param_info):
    plt.figure(figsize=(15, 10))
    plt.grid(visible = True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('Fraction of sites are methylated for that cell',fontsize=20)
    plt.ylabel('Extinction probability starting with a single cell',fontsize=20)

    for i, data_for_specific_M in enumerate(data):
        M = len(data_for_specific_M) - 1
        xs = [i / M for i in range(M + 1)]
        if i == len(data) - 1:
            label = 'limiting distribution'
        else:
            label = f'M = {M}'
        plt.plot(xs, data_for_specific_M, label = label)


    plt.text(0.21, 0.01, f'{param_info}')


    plt.legend(loc = 'lower left')

    plt.show()

