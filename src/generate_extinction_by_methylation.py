import constants
import models
import matplotlib.pyplot as plt


def plot_varying_methylation(data, scales, param_info = None):
    plt.figure(figsize=(12, 8))
    plt.grid(visible = True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('Fraction of sites are methylated for that cell',fontsize=20)
    plt.ylabel('Extinction probability starting with a single cell',fontsize=20)

    for i, data_for_methylation in enumerate(data):
        M = len(data_for_methylation) - 1
        xs = [i / M for i in range(M + 1)]
        label = f'r_mu = {scales[i]}'
        plt.plot(xs, data_for_methylation, label = label)


    plt.text(0.21, 0.01, f'{param_info}')


    plt.legend(loc = 'lower left')
    plt.show()




if __name__ == '__main__':
    base_params = constants.TOY_LIVING_PARAMS


    ratio = 2
    scales = [0.01, 0.3, 0.1, 1, 10]
    data = []
    r_mus = [scale for scale in scales]
    r_ums = [ratio * scale for scale in scales]
    for scale in scales:
        model_params = base_params[:]
        model_params[4] = ratio * scale
        model_params[5] = scale
        model = models.NoncollaborativeStochastic(model_params)
        extinction_rates = model.find_limit_extinction_probabilities()
        data.append(extinction_rates)
    
    param_info = ('Other parameters:\n' + 
                    'M: infinite\n' +
                    f'b: {base_params[1]}\nd_0: {base_params[2]}\n' + 
                    f'd_M: {base_params[3]}\nr_um / r_mu: {ratio}\n' + 
                    f'p: {base_params[6]}')
    

    plot_varying_methylation(data, scales, param_info=param_info)
 
