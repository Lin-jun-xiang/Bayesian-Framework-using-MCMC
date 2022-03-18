import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import Voronoi_tessellation as vt

prior_mass_flux = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\Prior_qc_v2_test.xlsx')['mass_flux']

posterior_v1 = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\Posterior_qc_v2.xlsx')['mass_flux']

posterior_v1_2 = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\Posterior_qc_v2-2.xlsx')['mass_flux']

# true_v1_file = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\\Obs_conc_v1.xlsx')
# vor_area = vt.voronoi(true_v1_file)
# c = true_v1_file['MINIT']
# true_v1 = 0
# for i, node in enumerate(vor_area):
#     true_v1 += 0.391245*c[i]*365*1e-3*vor_area[int(node)]

# true_v1 = true_v1/1000

# true_v1_2_file = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\\Obs_conc_v1-2.xlsx')
# vor_area = vt.voronoi(true_v1_2_file)
# c = true_v1_2_file['MINIT']
# true_v1_2 = 0
# for i, node in enumerate(vor_area):
#     true_v1_2 += 0.391245*c[i]*365*1e-3*vor_area[int(node)]

# true_v1_2 = true_v1_2/1000


sns.set()
plt.axvline(0.0105, color='red')
# plt.axvline(0.00966, color='blue')
sns.distplot(prior_mass_flux, label='prior', color='dodgerblue', hist=False, kde_kws={'linewidth':4})
sns.distplot(posterior_v1, label='posterior', color='coral', hist=False, kde_kws={'linewidth':4})
# sns.distplot(posterior_v1_2, label='15 sampling points', color='mediumorchid', hist=False, kde_kws={'linewidth':4})

plt.legend()
# plt.title('Case1-Mass flux')
plt.xlabel('mass flux(kg/(yr${⋅𝑚^2}$))', fontsize=15)
plt.show()
