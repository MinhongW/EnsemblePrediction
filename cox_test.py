import math
import numpy as np


# cox models
class CoxModel(object):
    # surv_baseline: baseline survival
    # beta_list: parameter estimates
    # x_mean_list: mean value of the variables
    # x_list: variables of the patient for risk calculation
    def __init__(self, surv_baseline, beta_list, x_mean_list, x_list):
        self.surv_baseline = surv_baseline
        self.beta = beta_list
        self.x = x_list
        self.x_mean = x_mean_list

    def cox_risk_eq1(self): # using equation 1
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
        return risk

# test swedish NDR
# 5 years survival baseline
surv_ndr = 0.97136
# parameter estimates:
# Diabetes duration
# Onset age of diabetes
# Log of TC:HDL
# Log of HbA1c(DCCT)
# Log of systolic BP
# Smoker
# Macroalbuminuria
# Previous CVD
beta_list_ndr = [0.08426,
                 0.04742,
                 0.80050,
                 1.27275,
                 1.20050,
                 0.56688,
                 0.41995,
                 1.25506]

# mean values of the variables listed above
x_mean_list_ndr = [28.014,
                   16.601,
                   1.1470,
                   2.0605,
                   4.8598,
                   0.1483,
                   0.1237,
                   0.0612]

# variables of an example patient
x_list_ndr_test = [30,
                   18,
                   np.log(4.55),
                   np.log(8),
                   np.log(150),
                   0,
                   1,
                   0]

# test the NDR model with the example patient
test_ndr = CoxModel(surv_ndr, beta_list_ndr, x_mean_list_ndr, x_list_ndr_test)
risk_ndr = test_ndr.cox_risk_eq1()
print(risk_ndr)

# test Fremantle Diabetes Study risk score
# 5 years survival baseline
surv_fds = 0.904

# parameters estimate:
# Age
# Male sex
# Prior CVD
# Ln (urinary albumin:creatinine ratio)
# Ln[HbA1c(%)]
# Ln(serum HDL-cholesterol)
# Southern European
# Indigenous Australian
beta_list_fds = [0.080, 0.335, 0.693, 0.872, 0.196, -0.608, -0.771, 1.269]
# mean values of the listed parameters
x_mean_list_fds = [64.069, 0.492, 0.323, 2.024, 1.135, 0.012, 0.183, 0.014]
# example 59 year-old patient
x_list_fds_test = [59, 1, 1, np.log(8.0), np.log(0.92), np.log(0.79), 1, 0]

# test the Fremantle model with the example patient
test_fds = CoxModel(surv_fds, beta_list_fds, x_mean_list_fds, x_list_fds_test)
risk_fds = test_fds.cox_risk_eq1()
print(risk_fds)

# test ADVANCE --> need double check the model
# surv_advance = 0.951044
# beta_list_advance = [0.06187, -0.4736, 0.08263, 0.00665, 0.38248, 0.60106, 0.09945, 0.19341, 0.12619, 0.24219]
# x_list_advance_test = [50, 1, 10, 100, 0, 0, 7.5, 2.05, 4.0, 0]
#
# test_advance = CoxModel(surv_advance, beta_list_advance, [], x_list_advance_test)
# risk_advance = test_advance.cox_risk_eq2()
# print(risk_advance)

# NZ DCS
