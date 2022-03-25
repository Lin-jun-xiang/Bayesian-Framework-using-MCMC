import pandas as pd
import xml.etree.ElementTree as ET
import ifm
import sys
from scipy import stats
import numpy as np

def get_fem_file(file):
    sys.path.append("C:\\Program Files\\DHI\\2020\\FEFLOW 7.3\\bin64")
    global doc
    doc = ifm.loadDocument(file)

def get_well_index(xml_file):
    """
    get well sampling data index including K, head, qc from xml file
    return index of elements "or" nodes
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    obs_index = []
    for child in root:
        obs_index.append(child.attrib['end'])

    return obs_index

def get_obs_data(obs_data_filepath):
    """
    return "true" data on sampling(observation)
    """
    obs_data = {}

    well_data = pd.read_excel(obs_data_filepath)
    for i in range(len(well_data)):
        obs_data[well_data['Node'][i]] = well_data['MINIT'][i]

    return obs_data

def get_sim_data(nodes, area, control_plane_area=1000):
    """
    Darcy flux - [m/d] to [m/yr] ; Concentration - [mg/l] to [kg/m^3] ; Area - [m^2]
    Mass discharge = [(Darcy flux)*365] * [(Concentration)*1e-3] * [(Area)] which unit is [kg/yr]
    return "simulated" data (head, mass concentration, darcy flux, mass discharge)
    """
    head, concentration, darcy_flux = {}, {}, {}
    mass_discharge = {}

    initialize_results()

    doc.startSimulator()
    for node in nodes:
        head[int(node)] = doc.getResultsFlowHeadValue(int(node)-1)
        darcy_flux[int(node)] = doc.getResultsXVelocityValue(int(node)-1)
        concentration[int(node)] = doc.getResultsTransportMassValue(int(node)-1)
        mass_discharge[int(node)] = darcy_flux[int(node)] * concentration[int(node)] * area[int(node)] * 365 * 1e-3

    mass_discharge['mass_discharge'] = sum([v for v in mass_discharge.values()])
    mass_discharge['mass_flux'] = mass_discharge['mass_discharge']/control_plane_area
    doc.stopSimulator()

    return {'head':head, 'concentration':concentration, 'darcy_flux':darcy_flux, 'mass_discharge':mass_discharge}

def set_ifm_K(zone_elements, theta_K):
    """
    Set conductivity(X=Y=10*Z) in each element; units - "m/d"
    Input conductivity units is "m/s", need to convert to "m/d"
    """
    if type(theta_K) == float:
        theta_K = theta_K * 24 * 3600

        for element in zone_elements:
            doc.setMatXConductivityValue3D(int(element)-1, theta_K)
            doc.setMatYConductivityValue3D(int(element)-1, theta_K)
            doc.setMatZConductivityValue3D(int(element)-1, theta_K*0.1)
    else:
        for i in range(len(zone_elements)):
            doc.setMatXConductivityValue3D(int(zone_elements[i])-1, theta_K[i] * 24 * 3600)
            doc.setMatYConductivityValue3D(int(zone_elements[i])-1, theta_K[i] * 24 * 3600)
            doc.setMatZConductivityValue3D(int(zone_elements[i])-1, theta_K[i] * 24 * 3600 * 0.1)

def initialize_results():
    """
    For transport simulation, need to initialize 'concentration' for avoid concentration accumulation.
    """
    [doc.setResultsTransportMassValue(node, 0) for node in range(doc.getNumberOfNodes())]

def covariance_matrix(n_obs, r=0.8, sigma=0.01):
    """
    Hyper-parameters: r, sigma. (reference: Bodin et al.,2012； Zheng et al.,2016)
    return the covariance matrix of data noise.
    """
    covMatrix = np.array([[1.0]*n_obs]*n_obs) # initialize

    for i in range(len(covMatrix)):
        for j in range(len(covMatrix)):
            covMatrix[i][j] = r**(abs(i-j))

    return sigma*covMatrix

def likelihood_calculate(obs_data, nodes, K_zone, theta, covMatrix, area, target='concentration'):
    """
    P(D|theta) = 1/((2π)^(N/2)*|C|^(1/2))*exp⁡((-1/2)*(d-f(θ))^T C^(-1) (d-f(θ))). (reference: Bodin et al.,2012； Zheng et al.,2016)
    left = 1/((2π)^(N/2)*|C|^(1/2))
    right = exp⁡((-1/2)*(d-f(θ))^T C^(-1) (d-f(θ)))
    return-1 P(D|theta) which is a value
    return-2 dictionary of simulate data(ex:concentration)
    """
    for i in range(len(theta)):
        set_ifm_K(K_zone[i], theta[i])

    sim_data = get_sim_data(nodes, area)

    residual = []
    for key in obs_data:
        residual.append(obs_data[key] - sim_data[target][key])

    likelihood_left =  1/((2*3.14159)**len(obs_data)*np.linalg.det(covMatrix))**0.5
    likelihood_right = np.math.exp(-0.5*np.array(residual).dot(np.linalg.inv(covMatrix)).dot(np.array(residual)))

    return {'likelihood':likelihood_left*likelihood_right, 'sim_data':sim_data}

def k_target_distribution(theta, mu=1e-3, sigma=1e-3):
    """
    prior conductivity distribution ~ N(mu, sigma)
    return a value which belong to the distr
    """
    return joint_distribution([stats.norm.pdf(np.math.log(theta_i), loc=np.math.log(mu**2/(mu**2+sigma**2)**0.5), scale=(np.math.log((sigma/mu)**2+1))**0.5) for theta_i in theta])

def proposal_distribution(theta, s=1):
    # return abs(stats.norm.rvs(loc=theta, scale=1e-4, size=s))
    return np.random.uniform(1e-6, 1e-3, size=2)

def proposal_calculate(target, condition):
    """
    Parameters
    ------
    target : theta_cur
    condition : theta_star

    """
    return joint_distribution(stats.norm.pdf(target, condition))

def joint_distribution(dist):
    """
    ln(results) = ln(P(x1)) + ln(P(x2)) +...
    return a value of joint distribution : P(x1)*P(x2)*...P(xn)
    """
    res = 0
    for v in dist:
        res += np.math.log(v)
    return np.math.exp(res)
