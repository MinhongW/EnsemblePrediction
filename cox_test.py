import math
import numpy as np


# cox models
class CoxModel(object):
    # surv_baseline: baseline survival
    # beta_list: parameter estimates
    # x_mean_list: mean value of the variables (used in eq1)
    # x_mean_adjusted: one parameter from the x mean list (used in eq2)
    # x_list: variables of the patient for risk calculation
    def __init__(self, surv_baseline, beta_list, x_mean_list, x_mean_adjusted, x_list):
        self.surv_baseline = surv_baseline
        self.beta = beta_list
        self.x = x_list
        self.x_mean = x_mean_list
        self.x_mean_adjusted = x_mean_adjusted

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
        a = a - self.x_mean_adjusted
        risk = 1 - self.surv_baseline ** math.exp(a)
        return risk

####################
#### swedish NDR ###
####################
# using equation 1 (cox_risk_eq1)
# ---- parameters are different from Steph's code, but I got parameters from the paper...
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
test_ndr = CoxModel(surv_ndr, beta_list_ndr, x_mean_list_ndr, [], x_list_ndr_test)
risk_ndr = test_ndr.cox_risk_eq1()
print(risk_ndr)

#################
### Fremantle ###
#################

# test Fremantle Diabetes Study risk score (same as steph's code)
# 5 years survival baseline
surv_fds = 0.904

# parameters estimate:
# Age
# Male sex
# Prior CVD
# Ln (urinary albumin:creatinine ratio)
# Ln[HbA1c(%)]
# Ln(serum HDL-cholesterol)
# if Southern European
# if Indigenous Australian
beta_list_fds = [0.080, 0.335, 0.693, 0.872, 0.196, -0.608, -0.771, 1.269]
# parameter from mean value
x_mean_adjusted_fds = 7.371
# example 59 year-old patient
x_list_fds_test = [59, 1, 1, np.log(8.0), np.log(0.92), np.log(0.79), 1, 0]

# test the Fremantle model with the example patient
test_fds = CoxModel(surv_fds, beta_list_fds, [], x_mean_adjusted_fds, x_list_fds_test)
risk_fds = test_fds.cox_risk_eq2()
print(risk_fds)

###############
### ADVANCE ###
###############
# test ADVANCE (from Steph's code)
# this is the surv baseline for 5-year risk, 4-year risk is also provided in the paper
surv_advance = 0.938805

# parameters:
# age at diagnosis (per 1 year increase)
# known duration of diabetes (per 1 year increase)
# pulse pressure (per 1 mmhg increase)
# retinopathy (yes or no)
# atrial fibrillation (old/present vs absent)
# hba1c (per 1% increase)
# log of urinary albumin/creatinine ratio (per 1 log mg/g increase)
# non-hdl cholesterol (per 1 mmol/l increase)
# treated hypertension (yes or no)
# sex (women vs men)
beta_list_advance = [0.06187, 0.08263, 0.00665, 0.38248, 0.60106, 0.09945, 0.19341, 0.12619, 0.24219, -0.4736]
x_mean_adjusted_advance = 6.55666

x_list_advance_test = [50, 10, 100, 0, 0, 8, 2.05, 1, 1, 0]
test_advance = CoxModel(surv_advance, beta_list_advance, [], x_mean_adjusted_advance, x_list_advance_test)
risk_advance = test_advance.cox_risk_eq2()
print(risk_advance)

# NZ DCS
surv_nzdcs = 0.8156

beta_list_nzdcs = [0.04093, -0.16542, 0.10296, 0.23131, 0.05874, 0.00290, -0.09536, 0.22218, 0.02015, 0.21857, 0.69659, 0.05795]
x_mean_adjusted_nzdcs = 3.62948747

x_list_nzdcs_test = [50, 0, 7, 7, 8, 0, 0, 1, 8, 0, 0, 10]

test_nzdcs = CoxModel(surv_nzdcs, beta_list_nzdcs, [], x_mean_adjusted_nzdcs, x_list_nzdcs_test)
risk_nzdcs = test_nzdcs.cox_risk_eq2()
print(risk_nzdcs)

###########
### CHS ###
###########
# CHS includes two separated models for male/female
# using eq1 function here
# for female
surv_chs_f = 0.8653247
beta_list_chs_f = [0.0526546, 0.257817, 0.49725, 0.01436, 0.041323, -0.0620634, 0.358296, 0.5346682]
x_mean_list_chs_f = [72.6, 0, 0, 141/10, 5.69, 1.35, 0, 0]

x_list_chs_f_test = [60, 0, 0, 150/10, 6, 1, 0, 0]

test_chs_f = CoxModel(surv_chs_f, beta_list_chs_f, x_mean_list_chs_f, [], x_list_chs_f_test)
risk_chs_f = test_chs_f.cox_risk_eq1()
print(risk_chs_f)

# for male
surv_chs_m = 0.7709767
beta_list_chs_m = [0.0526546, 0.257817, 0.49725, 0.01436, 0.041323, 0.0620634, 0.358296, 0.5346682]
x_mean_list_chs_m = [73, 0, 0, 14, 5.05, 1.15, 0, 0]

x_list_chs_m_test = [60, 0, 0, 150/10, 6, 1, 0, 0]

test_chs_m = CoxModel(surv_chs_m, beta_list_chs_m, x_mean_list_chs_m, [], x_list_chs_m_test)
risk_chs_m = test_chs_m.cox_risk_eq1()
print(risk_chs_m)
