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








#def processOCV(data,cellID,minV,maxV,savePlots)
