import scipy.io, math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.linalg as la
import scipy.optimize
from scipy.optimize import nnls,brentq
from scipy.interpolate import interp1d
from joblib import dump,load
from processOCV import processOCV

from oct2py import octave
octave.addpath('~/github/bms-specialization/course-2/scripts/')

# --------------------------------------------------------------------
# function processDynamic
#
# Technical note: PROCESSDYNAMIC assumes that specific Arbin test
# scripts have been executed to generate the input files.
# "makeMATfiles.m" converts the raw Excel data files into "MAT" format
# where the MAT files have fields for time, step, current, voltage,
# chgAh, and disAh for each script run.
#
# The results from three scripts are required at every temperature.
# The steps in each script file are assumed to be:
#   Script 1 (thermal chamber set to test temperature):
#     Step 1: Rest @ 100# SOC to acclimatize to test temperature
#     Step 2: Discharge @ 1C to reach ca. 90# SOC
#     Step 3: Repeatedly execute dynamic profiles (and possibly
#             intermediate rests) until SOC is around 10#
#   Script 2 (thermal chamber set to 25 degC):
#     Step 1: Rest ca. 10# SOC to acclimatize to 25 degC
#     Step 2: Discharge to min voltage (ca. C/3)
#     Step 3: Rest
#     Step 4: Constant voltage at vmin until current small (ca. C/30)
#     Steps 5-7: Dither around vmin
#     Step 8: Rest
#   Script 3 (thermal chamber set to 25 degC):
#     Step 2: Charge @ 1C to max voltage
#     Step 3: Rest
#     Step 4: Constant voltage at vmax until current small (ca. C/30)
#     Steps 5-7: Dither around vmax
#     Step 8: Rest
#
# All other steps (if present) are ignored by PROCESSDYNAMIC. The time
# step between data samples must be uniform -- we assume a 1s sample
# period in this code
#
# The inputs:
# - data: An array, with one entry per temperature to be processed.
#         One of the array entries must be at 25 degC. The fields of
#         "data" are: temp (the test temperature), script1,
#         script 2, and script 3, where the latter comprise data
#         collected from each script.  The sub-fields of these script
#         structures that are used by PROCESSDYNAMIC are the vectors:
#         current, voltage, chgAh, and disAh
# - model: The output from processOCV, comprising the OCV model
# - numpoles: The number of R-C pairs in the model
# - doHyst: 0 if no hysteresis model desired; 1 if hysteresis desired
#
# The output:
# - model: A modified model, which now contains the dynamic fields
#         filled in.
  # used by fminbnd later on
  # options=optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEval',100000, ...
  #   'MaxIter',1e6,'Jacobian','Off'); # for later optimization
  # options=optimset('TolX',0.1,'TolFun',1e-2,'MaxFunEval',40, ...
  #   'MaxIter',20,'Jacobian','Off'); # for later optimization


  # function model = processDynamic(data,model,numpoles,doHyst)
  #   global bestcost


class ProcessDynamic:

    def __init__(self, \
            data_dir='P14_DYN/', \
            model_dir='data/', \
            num_poles=1, \
            do_hyst=1):
        self.model = scipy.io.loadmat('../data/p14model-ocv-mat7.mat', simplify_cells=True)['model']
        self.data_dir = data_dir
        self.model_dir = model_dir
        self.temps = [5,25,45]
        self.magnitudes = [30,50,50]
        self.scripts = ['script1','script2','script3']
        self.num_poles = num_poles
        self.do_hyst = do_hyst
        self.data = dict()
        self.bestcost = np.inf

    def load_data(self,verbose=False):

        for temp,mag in zip(self.temps,self.magnitudes):
            if temp < 0:
              filename = self.model_dir + self.data_dir + 'P14_DYN_%d_N%d.mat' % (mag,abs(temp))
            else:
              filename = self.model_dir + self.data_dir + 'P14_DYN_%d_P%d.mat' % (mag,abs(temp))

            self.data[temp] = scipy.io.loadmat(filename, simplify_cells=True)['DYNData']
            if verbose:
              print('Loaded: ',filename)

        for key in self.data.keys():
            self.model[key] = dict()
            self.model[key]['GParam'] = 0
            self.model[key]['RCParam'] = [0]

    def process_DYN_step1(self):

        # Step 1: Compute capacity and coulombic efficiency for every test
        self.temps = list(self.data.keys())

        if 25 not in self.temps:
            raise FileNotFoundError

        self.temps.remove(25)
        self.temps = [25] + self.temps

        for temp in self.temps:
            total_discharge = 0
            total_charge = 0

            for script in self.scripts:
                total_discharge += self.data[temp][script]['disAh'][-1]
                total_charge += self.data[temp][script]['chgAh'][-1]

            eta_t = total_discharge / total_charge

            if temp == 25:
              self.eta_25 = eta_t

            Q_t = self.data[temp]['script1']['disAh'][-1] \
            + self.data[temp]['script2']['disAh'][-1] \
            - self.eta_25*self.data[temp]['script1']['chgAh'][-1] \
            - self.eta_25*self.data[temp]['script2']['chgAh'][-1]
            print(temp,Q_t,eta_t)

            self.model[temp]['Q'] = Q_t
            self.model[temp]['eta'] = eta_t

    def process_DYN_step2(self):

        # Step 2: Compute OCV for "discharge portion" of test
        model = processOCV(data_dir='../data/P14_OCV/')
        model.run()

        for temp in self.temps:

            eta_param = self.model[temp]['eta']
            charge_indices = np.where(self.data[temp]['script1']['current'] < 0)
            corrected_current = self.data[temp]['script1']['current']
            corrected_current[charge_indices] = eta_param*corrected_current[charge_indices]

            self.model[temp]['Z'] = 1 - np.cumsum(corrected_current)/int(self.model[temp]['Q']*3600)
            self.model[temp]['OCV'] = model.SOC_temp_to_OCV(self.model[temp]['Z'],temp)

    def process_DYN_step3(self,do_hyst):

        for temp in self.temps:
            print('Processing temperature: ',temp)

            if do_hyst:
                function = lambda x_G: self.optfn(x_G,temp,do_hyst)
                self.model[temp]['G'] = scipy.optimize.fminbound(function,1,250,xtol=0.1,maxfun=40)
            else:
                GParam = 0
                self.optfn(GParam,temp,do_hyst)

    def optfn(self,GParam,temp,do_hyst):

        self.model[temp]['G'] = GParam
        #print('Updated G value, optfn', GParam, 'temp',temp,'do_hyst',do_hyst)
        cost = self.minfn(temp,do_hyst,GParam)
        #print(GParam,cost,self.bestcost)

        if cost < self.bestcost:
            self.bestcost = cost
            print('The model created for this value of gamma is the best ESC model yet! GParam = ',GParam)

        return cost


    def minfn(self,temp,do_hyst,GParam):
        #
        # model = processOCV(data_dir='../data/P14_OCV/')
        # model.run()

        numfiles = len(self.temps)
        xplots = math.ceil(math.sqrt(numfiles))
        yplots = math.ceil(numfiles/xplots)
        rmserr = np.zeros((1,xplots*yplots))

        #G = self.model[temp]['GParam']
        G = GParam
        #print('Updated G value, minfn', G, 'temp',temp,'do_hyst',do_hyst)
        Q = self.model[temp]['Q']
        eta = self.model[temp]['eta']
        # RC = self.model[temp]['RCParam']
        # numpoles = len(RC)
        # # print(temp,G)
        # # print(temp,Q)
        # # print(temp,eta)
        # # print(temp,RC)
        # # print(temp,numpoles)

        current = self.data[temp]['script1']['current']
        voltage = self.data[temp]['script1']['voltage']
        #time = np.ar
        charge_indices = np.where(self.data[temp]['script1']['current'] < 0)
        corrected_current = self.data[temp]['script1']['current']
        corrected_current[charge_indices] = eta*corrected_current[charge_indices]

        h = 0*current
        s = 0*current

        fac = np.exp(-abs(G*corrected_current/(3600*Q)))
        #print('fac ',fac[1000], 'G ',G)

        for i in range(1,len(current)):
            h[i] = fac[i-1]*h[i-1] + (fac[i-1]-1)*np.sign(current[i-1])
            s[i] = np.sign(current[i])
            if abs(current[i]) < Q/100:
                s[i] = s[i-1]

        # First modeling step: Compute error with model = OCV only
        vest1 = self.model[temp]['OCV'][0]
        verr = voltage - vest1
        self.model[temp]['vest1'] = vest1
        self.model[temp]['verr1'] = verr
        self.model[temp]['current'] = current
        self.model[temp]['h'] = h
        self.model[temp]['s'] = s
        self.model[temp]['corrected_current'] = corrected_current

        # Second modeling step: Compute time constants in "A" matrix
        num_poles = self.num_poles

        while True:
            A = octave.SISOsubid(-np.diff(verr),np.diff(corrected_current),num_poles)
            eigA = la.eig(np.identity(1)*A)[0]
            eigA = eigA[eigA == np.conj(eigA)]
            eigA = eigA[eigA > 0 and eigA < 1]
            okpoles = len(eigA)
            num_poles = num_poles+1
            if okpoles >= self.num_poles:
                break

        RCfact = np.sort(eigA)
        # COME BACK placeholder for indexing
        RC = -1./np.log(RCfact)

        vrcRaw = np.zeros((self.num_poles,len(h)))

        self.model[temp]['vrcRaw'] = vrcRaw

        for i in range(1,len(current)):
            vrcRaw[:,i] = np.diag(RCfact)*vrcRaw[:,i-1] + (1-RCfact)*current[i-1]

        if do_hyst:
            #print('hyst on')
            H = [h,s,-current,-vrcRaw[0]]
            H = np.vstack((H[0],H[1],H[2],H[3]))
            W = nnls(H.T,verr)[0]
            M = W[0]
            M0 = W[1]
            R0 = W[2]
            Rfact = W[3]
            #print('look here',M,M0,R0,Rfact)
            #print(H[0][-1],H[1][-1],H[2][-1],H[3][-1])

        else:
            #print('hyst off')
            H = [-current,-vrcRaw.T]
            H = np.vstack((H[0],H[1].T[0]))
            W = nnls(H.T, verr)[0]
            M=0
            M0=0
            R0 = W[0]
            Rfact = W[1]

        vest2 = vest1 + M*h + M0*s - R0*current - vrcRaw[0]*Rfact
        verr = voltage - vest2

        # Set params
        self.model[temp]['R0'] = R0
        self.model[temp]['M0'] = M0
        self.model[temp]['M'] = M
        self.model[temp]['RC'] = RC
        self.model[temp]['Rfact'] = Rfact

        self.model[temp]['test1'] = R0*current
        self.model[temp]['test2'] = vrcRaw*Rfact
        self.model[temp]['vest2'] = vest2
        self.model[temp]['verr2'] = verr

        OCVstatic = processOCV(data_dir='../data/P14_OCV/')
        OCVstatic.run()

        v1 = OCVstatic.SOC_temp_to_OCV(0.95,temp)[0]
        v2 = OCVstatic.SOC_temp_to_OCV(0.05,temp)[0]
        N1 = np.where(voltage<v1)[0]
        N2 = np.where(voltage<v2)[0]

        if N1.size == 0:
            N1 = 1
        else:
            N1 = N1[0]

        if N2.size == 0:
            N2 = len(verr)
        else:
            N2 = N2[0]

        rmserr = np.sqrt(np.mean(verr[N1:N2]**2))
        cost = np.sum(rmserr)

        print('RMS error for present value of gamma = %f (mV)' % (np.round(cost*1000,2)))

        return cost
