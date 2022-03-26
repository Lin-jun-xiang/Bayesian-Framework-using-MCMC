# For Demo case (heterogeneous, random field)
import Mc_module as MC
import numpy as np
import pandas as pd
import time
import Voronoi_tessellation
import RandomField as RF

def getPositiondata(eleFile):
    """
    Get coordinates of "elements" from fem file, then
    return the pos for random field generator
    """
    pos = (eleFile["X"], eleFile["Y"], eleFile["Z"])

    return pos

def prior_ensemble():
    prior_gs_mean, prior_gs_var, prior_gs_ls = [], [], []
    realizations = []
    prior = []
    for i in range(nRealizations):
        # Only take rd_mean in prior
        prior_gs_mean.append(rfg.gs_meanDistr()[0])
        prior_gs_var.append(rfg.gs_varDistr())
        prior_gs_ls.append(rfg.gs_lenScaleDistr())
        realizations.append(rfg.uncondFieldGenerator(pos,
                                                     prior_gs_mean[i],
                                                     prior_gs_var[i],
                                                     prior_gs_ls[i]))

        MC.set_ifm_K(K_zone, realizations[i])
        prior.append(MC.get_sim_data(obs_nodes, area=voronoi_area)[output_target])
        print('prior-', i+1)

    return prior

def create_markov_chain(n):
    chain = []
    # Initialize the uncertainty parameters by sampling from proposal distr
    # Notice the gs_var and gs_lenScale which itself is uniform distr (eg proposal distr)
    chain.append([MC.gs_mean_proposal_distribution(),
                 rfg.gs_varDistr(),
                 rfg.gs_lenScaleDistr()])

    rejection_rate = 0

    covMatrix = MC.covariance_matrix(n_obs=len(obs_data))

    posterior = []
    for t in range(n-1):
        theta_cur = chain[-1]
        realization_cur = rfg.uncondFieldGenerator(pos,
                                                   theta_cur[0],
                                                   theta_cur[1],
                                                   theta_cur[2])
        likelihood_cur = MC.likelihood_calculate(obs_data, obs_nodes, K_zone, realization_cur, covMatrix, area=voronoi_area)

        theta_star = [MC.gs_mean_proposal_distribution(),
                      rfg.gs_varDistr(),
                      rfg.gs_lenScaleDistr()]
        realization_star = rfg.uncondFieldGenerator(pos,
                                                    theta_star[0],
                                                    theta_star[1],
                                                    theta_star[2])
        likelihood_star = MC.likelihood_calculate(obs_data, obs_nodes, K_zone, realization_star, covMatrix, area=voronoi_area)

        u = np.random.uniform(size=1)[0]

        acceptance_rate = min(1, likelihood_star['likelihood']*MC.gs_parameters_target_distribution(theta_star)*MC.proposal_calculate(theta_cur, theta_star)/\
            MC.gs_parameters_target_distribution(theta_cur)/MC.proposal_calculate(theta_star, theta_cur)/likelihood_cur['likelihood'])
        # acceptance_rate = min(1, likelihood_star['likelihood']/likelihood_cur['likelihood'])

        if u <= acceptance_rate:
            chain.append(theta_star)
            posterior.append(likelihood_star['sim_data'][output_target])
        else:
            chain.append(theta_cur)
            posterior.append(likelihood_cur['sim_data'][output_target])
            # If rejection_rate too big, can adjust the "u"
            rejection_rate += 1
        print('chain-', t+1)

    return {'chain':chain, 'posterior':posterior, 'rejection':rejection_rate/(n-1)}

if __name__ == "__main__":
    time_start = time.time()

    MC.get_fem_file('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\True_transport_RandomField.fem')

    K_zone = [e+1 for e in range(MC.doc.getNumberOfElements())]

    eleFile = pd.read_excel("C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\True_transport_RandomField_Kxx.fem.xlsx")
    pos = getPositiondata(eleFile)

    obs_data_file = 'C:\\Users\\JunXiang\Desktop\\傑明工程\\fem\\excel\\Obs_conc_RandomField.xlsx'
    obs_data = MC.get_obs_data(obs_data_file)

    obs_nodes = list(obs_data.keys())

    voronoi_area = Voronoi_tessellation.voronoi(obs_data=pd.read_excel(obs_data_file))
    control_plane_area = 1000

    output_target = 'mass_discharge'

    nRealizations = 100

    rfg = RF.RandomFieldGenerator()

    prior = prior_ensemble()

    markov_chain = create_markov_chain(n=nRealizations)
    burn_in_period = 0

    posterior = markov_chain['posterior'][burn_in_period-1:]

    time_end = time.time()

    print('rejection rate=', markov_chain['rejection'])
    print('time=', time_end-time_start)
    # pd.DataFrame(markov_chain['chain']).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Theta_v1.xlsx')
    pd.DataFrame(prior).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Prior_qc_randomField.xlsx')
    pd.DataFrame(posterior).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Posterior_qc_randomField.xlsx')
    
    import winsound
    winsound.PlaySound('SystemHand', winsound.SND_ALIAS)
