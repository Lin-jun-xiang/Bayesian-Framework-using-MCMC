# For Minshuan case (heterogeneous, random field)
import Mc_module as MC
from scipy import stats
import numpy as np
import pandas as pd
import time
import Voronoi_tessellation
import RandomField as RF

def getPositiondata():
    pos = ([], [], [])
    for i in range(MC.doc.getNumberOfNodes()):
        pos[0].append(MC.doc.getX(i))
        pos[1].append(MC.doc.getY(i))
        pos[2].append(MC.doc.getZ(i))

    return pos

def prior_ensemble(mu=1e-3, sigma=1e-3):
    rfg = RF.RandomFieldGenerator()
    prior_gs_mean, prior_gs_var, prior_gs_ls = [], [], []
    realizations = []
    for i in range(nRealizations):
        # Only take rd_mean in prior
        prior_gs_mean.append(rfg.gs_meanDistr()[0])
        prior_gs_var.append(rfg.gs_varDistr())
        prior_gs_ls.append(rfg.gs_lenScaleDistr())
        realizations.append(rfg.uncondFieldGenerator(pos,
                                                     prior_gs_mean[i],
                                                     prior_gs_var[i],
                                                     prior_gs_ls[i])
    prior = []
    for i in range(len(prior)):
        MC.set_ifm_K(K1_zone_elements, np.math.exp(prior_theta_K1[i]))
        MC.set_ifm_K(K2_zone_elements, np.math.exp(prior_theta_K2[i]))
        prior.append(MC.get_sim_data(nodes, area=voronoi_area)[output_target])
        print('prior-', i+1)

    return prior

def create_markov_chain(n):
    chain = []
    chain.append(MC.proposal_distribution(theta=1e-3, s=2))
    posterior = []
    rejection_rate = 0

    covMatrix = MC.covariance_matrix(n_obs=len(obs_data))

    for t in range(n-1):
        theta_cur = chain[-1]

        likelihood_cur = MC.likelihood_calculate(obs_data, nodes, K_zone, theta_cur, covMatrix, area=voronoi_area)

        theta_star = MC.proposal_distribution(theta=theta_cur, s=2)

        likelihood_star = MC.likelihood_calculate(obs_data, nodes, K_zone, theta_star, covMatrix, area=voronoi_area)

        u = np.random.uniform(size=1)[0]

        acceptance_rate = min(1, likelihood_star['likelihood']*MC.k_target_distribution(theta_star)*MC.proposal_calculate(theta_cur, theta_star)/\
            MC.k_target_distribution(theta_cur)/MC.proposal_calculate(theta_star, theta_cur)/likelihood_cur['likelihood'])
        # acceptance_rate = min(1, likelihood_star['likelihood']/likelihood_cur['likelihood'])

        if u <= acceptance_rate:
            chain.append(theta_star)
            posterior.append(likelihood_star['sim_data'][output_target])
        else:
            chain.append(theta_cur)
            posterior.append(likelihood_cur['sim_data'][output_target])
            rejection_rate += 1
        print('chain-', t+1)

    return {'chain':chain, 'posterior':posterior, 'rejection':rejection_rate/(n-1)}

if __name__ == "__main__":
    time_start = time.time()

    MC.get_fem_file('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\True_transport.fem')

    K1_zone_elements = pd.read_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\K1_ele.xlsx')['Element']
    K2_zone_elements = pd.read_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\K2_ele.xlsx')['Element']
    K_zone = [K1_zone_elements, K2_zone_elements]

    nodes = MC.get_well_index("C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\xml\\well_sampling2.xml")
    # elements = MC.get_well_index("C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\xml\\well_K.xml")
    obs_data_file = 'C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Obs_conc_v1-2.xlsx'
    obs_data = MC.get_obs_data(obs_data_file)

    voronoi_area = Voronoi_tessellation.voronoi(obs_data=pd.read_excel(obs_data_file))
    control_plane_area = 1000

    output_target = 'mass_discharge'

    nRealizations = 100

    prior = prior_ensemble()

    markov_chain = create_markov_chain(nRealizations)
    burn_in_period = 0

    posterior = markov_chain['posterior'][burn_in_period-1:]

    time_end = time.time()

    print('rejection rate=', markov_chain['rejection'])
    print('time=', time_end-time_start)
    # pd.DataFrame(markov_chain['chain']).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Theta_v1.xlsx')
    pd.DataFrame(prior).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Prior_qc_v1-2.xlsx')
    pd.DataFrame(posterior).to_excel('C:\\Users\\JunXiang\\Desktop\\傑明工程\\fem\\excel\\Posterior_qc_v1-2.xlsx')
    
    import winsound
    winsound.PlaySound('SystemHand', winsound.SND_ALIAS)