import math
import numpy as np


# cox models
class CoxModel(object):
    def __init__(self, surv_baseline, beta_list, x_list, x_mean_list):
        self.surv_baseline = surv_baseline
        self.beta = beta_list
        self.x = x_list
        self.x_mean = x_mean_list

    def cox_risk(self): # using equation 1
        a = 0
        for beta, x, x_mean in zip(self.beta, self.x, self.x_mean):
            a = a + beta * (x - x_mean)
        risk = 1 - self.surv_baseline ** math.exp(a)
        return risk

    def cox_risk_eq2(self): # using equation 2
        a = 0
        for beta, x in zip(self.beta, self.x):
            a = a + beta * x
        risk = 1 - self.surv_baseline ** math.exp(a)

# test swedish NDR
surv5_ndr = 0.97136
beta_list_ndr = [0.08426,
                 0.04742,
                 0.80050,
                 1.27275,
                 1.20050,
                 0.56688,
                 0.41995,
                 1.25506]
x_mean_list_ndr = [28.014,
                   16.601,
                   1.1470,
                   2.0605,
                   4.8598,
                   0.1483,
                   0.1237,
                   0.0612]
x_list_ndr_test = [30,
                   18,
                   np.log(4.55),
                   np.log(8),
                   np.log(150),
                   0,
                   1,
                   0]
test_ndr = CoxModel(surv5_ndr,beta_list_ndr,x_list_ndr_test,x_mean_list_ndr)
risk_ndr = test_ndr.cox_risk()
print(risk_ndr)

# test Fremantle Diabetes Study risk score
surv_fds = 0.904
beta_list_fds = [0.080, 0.335, 0.693, 0.872, 0.196, -0.608, -0.771, 1.269]
x_mean_list_fds = [64.069, 0.492, 0.323, 2.024, 1.135, 0.012, 0.183, 0.014]
x_list_fds_test = [59, 1, 1, np.log(8.0), np.log(0.92), np.log(0.79), 1, 0]

test = CoxModel(surv_fds,beta_list_fds,x_list_fds_test,x_mean_list_fds)
risk_fds = test.cox_risk()
print(risk_fds)

# test ADVANCE


