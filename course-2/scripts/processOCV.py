import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# Technical note: PROCESSOCV assumes that specific Arbin test scripts
# have been executed to generate the input files. "makeMATfiles.m"
# converts the raw Excel data files into "MAT" format, where the MAT
# files have fields for time, step, current, voltage, chgAh, and disAh
# for each script run.
#
# The results from four scripts are required at every temperature.
# The steps in each script file are assumed to be:
#   Script 1 (thermal chamber set to test temperature):
#     Step 1: Rest @ 100% SOC to acclimatize to test temperature
#     Step 2: Discharge @ low rate (ca. C/30) to min voltage
#     Step 3: Rest ca. 0%
#   Script 2 (thermal chamber set to 25 degC):
#     Step 1: Rest ca. 0% SOC to acclimatize to 25 degC
#     Step 2: Discharge to min voltage (ca. C/3)
#     Step 3: Rest
#     Step 4: Constant voltage at vmin until current small (ca. C/30)
#     Steps 5-7: Dither around vmin
#     Step 8: Rest
#     Step 9: Constant voltage at vmin for 15 min
#     Step 10: Rest
#   Script 3 (thermal chamber set to test temperature):
#     Step 1: Rest at 0% SOC to acclimatize to test temperature
#     Step 2: Charge @ low rate (ca. C/30) to max voltage
#     Step 3: Rest
#   Script 4 (thermal chamber set to 25 degC):
#     Step 1: Rest ca. 100% SOC to acclimatize to 25 degC
#     Step 2: Charge to max voltage (ca. C/3)
#     Step 3: Rest
#     Step 4: Constant voltage at vmax until current small (ca. C/30)
#     Steps 5-7: Dither around vmax
#     Step 8: Rest
#     Step 9: Constant voltage at vmax for 15 min
#     Step 10: Rest
#
# All other steps (if present) are ignored by PROCESSOCV.  The time
# step between data samples is not critical since the Arbin integrates
# ampere-hours to produce the two Ah columns, and this is what is
# necessary to generate the OCV curves.  The rest steps must
# contain at least one data point each.


class processOCV:

    def __init__(self,data_dir='data/P14_OCV/'):
        self.data_dir = data_dir
        self.data = dict()
        self.temps = [-25,-15,-5,5,15,25,35,45]
        self.scripts = ['script1','script2','script3','script4']
        self.SOC = np.arange(0,1,0.005)
        self.model = dict()

    def load_data(self,verbose=False):

        for temp in self.temps:
            if temp < 0:
                filename = self.data_dir + 'P14_OCV_N%d.mat' % (abs(temp))
            else:
                filename = self.data_dir + 'P14_OCV_P%d.mat' % (abs(temp))

            self.data[temp] = scipy.io.loadmat(filename, simplify_cells=True)['OCVData']
            if verbose:
                print('Loaded: ',filename)

        for key in self.data.keys():
            self.model[key] = dict()

    def process_OCV(self):

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

            df_1 = pd.DataFrame(self.data[temp]['script1'])

            # Beginning of discharge
            IR1Da = df_1[df_1['step'] == 1].tail(1).voltage.values[0] \
                    - df_1[df_1['step'] == 2].head(1).voltage.values[0]

            # End of discharge
            IR2Da = df_1[df_1['step'] == 3].head(1).voltage.values[0] \
                    - df_1[df_1['step'] == 2].tail(1).voltage.values[0]

            IndC = np.where(self.data[temp]['script3']['step'] == 2)[0]

            # Beginning of charge
            IR1Ca = self.data[temp]['script3']['voltage'][IndC[0]] \
                - self.data[temp]['script3']['voltage'][IndC[0]-1]

            # End of charge
            IR2Ca = self.data[temp]['script3']['voltage'][IndC[-1]] \
                - self.data[temp]['script3']['voltage'][IndC[-1]+1]

            # These are just some sanity checks to keep these correction factors
            # from getting too large in case of data issues
            IR1D = min(IR1Da,2*IR2Ca)
            IR2D = min(IR2Da,2*IR1Ca)
            IR1C = min(IR1Ca,2*IR2Da)
            IR2C = min(IR2Ca,2*IR1Da)

            # Discharge voltage & SOC curves
            IndD = np.where(self.data[temp]['script1']['step'] == 2)[0]
            blend = np.linspace(0,1,len(IndD))
            IRblend = IR1D + (IR2D-IR1D)*blend
            disV = self.data[temp]['script1']['voltage'][IndD] + IRblend

            # Z = State of Charge
            disZ = 1 - self.data[temp]['script1']['disAh'][IndD]/Q_t
            disZ = disZ + (1 - disZ[0])

            assert disZ.shape == disV.shape

            # Charge voltage & SOC curves
            blend = np.linspace(0,1,len(IndC))
            IRblend = IR1C + (IR2C-IR1C)*blend
            chgV = self.data[temp]['script3']['voltage'][IndC] - IRblend
            chgZ = self.data[temp]['script3']['chgAh'][IndC]/Q_t
            chgZ = chgZ - chgZ[0]

            assert chgZ.shape == chgV.shape

            charge_interpolate = interp1d(chgZ,chgV)
            discharge_interpolate = interp1d(disZ,disV)

            deltaV50 = charge_interpolate(0.5) - discharge_interpolate(0.5)

            # Select points on charge curve where SOC < 50%.
            ind = np.where(chgZ < 0.5)[0]
            vChg = chgV[ind] - chgZ[ind]*deltaV50
            zChg = chgZ[ind]

            # Select points on charge curve where SOC > 50%.
            ind = np.where(disZ > 0.5)[0]
            vDis = (disV[ind] + ((1 - disZ[ind])*deltaV50))[::-1]
            zDis = (disZ[ind])[::-1]

            rawocv_interpolation = interp1d(np.concatenate([zChg,zDis]), \
                                    np.concatenate([vChg,vDis]), \
                                    kind='linear',fill_value="extrapolate")
            rawocv = rawocv_interpolation(self.SOC)

            self.model[temp]['SOC'] = self.SOC
            self.model[temp]['rawocv'] = rawocv
            self.model[temp]['Q'] = Q_t
            self.model[temp]['eta'] = eta_t

    def create_lookup_table(self):
        # Use only data where T > 0C as cold temperate data
        temps = [x for x in self.data.keys() if x > 0]
        num_temps = len(temps)
        ocv_size = self.model[25]['rawocv'].shape[0]

        ocv_arr = np.ones((num_temps,ocv_size))
        temp_arr = np.ones((num_temps,1))

        for ind,temp in enumerate(temps):
            temp_arr[ind] = temp
            ocv_arr[ind] = self.model[temp]['rawocv']

        temp_arr = np.hstack((np.ones((len(temp_arr),1)),temp_arr))

        # Solve linear set of equations, A*x = y
        A = temp_arr
        y = ocv_arr

        self.OCV0, self.OCVrel = np.linalg.lstsq(A, y, rcond=None)[0]

    def run(self):

        self.load_data()
        self.process_OCV()
        self.create_lookup_table()

    def SOC_temp_to_OCV(self,z,temp,two_decimals=True):

        OCV_T = self.OCV0 + temp*self.OCVrel
        OCV_curve = interp1d(self.SOC, \
            OCV_T, \
            kind='linear', \
            fill_value='extrapolate')

        estimated_OCV = OCV_curve(z)

        if two_decimals == True:
            estimated_OCV = np.round(estimated_OCV,2)

        print('Estimated Open Circuit Voltage (T=%d C,SOC=%f): %f' % (temp,z,estimated_OCV))

        return estimated_OCV, OCV_T


if __name__ == '__main__':

    model = processOCV()
    model.run()

    OCV = model.SOC_temp_to_OCV(0.5,25)
