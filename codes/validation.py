from pmethod import *
from pkmethod import *

# Model validation for p-method using eigenvalue analysis (see Handout 1 for further informations)
# Method is validated by comparings results with Hodges & Pierces (Fig. 5.3 & 5.4 [p.186-187])
# -> This code has been validated and gives the same flutter speed and flutter frequency than the Hodges & Pierces book 

""" p_method = pmethod(a=-1/5, e=-1/20, mu=20, r=np.sqrt(6/25), sigma=2/5, V_vec=np.linspace(0.01, 3, 100))
p_method.run()
p_method.plot_results() """

# Model validation for p-k method using (see Handout 2 for further informations)
# Method is validated by comparings results with Hodges & Pierces (Fig. 5.19 & 5.20 [p.228-229])
# -> This code has been validated and gives the same flutter speed and flutter frequency than the Hodges & Pierces book 

pk_method = pkmethod(a=-1/5, e=-1/10, mu=20, r=np.sqrt(6/25), sigma=2/5, V_vec=np.linspace(0.01, 3, 500))
pk_method.run()
pk_method.plot_results()
pk_method.find_flutter()

