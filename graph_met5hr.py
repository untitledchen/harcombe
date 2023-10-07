import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# seed = random.randrange(1000)
# #
# lst = [('culture', 'met', 'rep', 'frac_slope')]
# mets = list(range(20, 40, 1))
# for met in mets:
#     run(seed, "mono", 10, (0.01, 0), 10, (met, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
#     #run(seed, "co", 5, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
#
#     run_calc_frac(f'hcb_sim_mono_{seed}_met{met}.csv', 1000, 5)
#
#     data1 = pd.read_csv(f'frac_hcb_sim_mono_{seed}_met{met}.csv', header=2, na_filter=False)
#
#     yvar = 'frac'
#     data1[f'log10_{yvar}'] = np.log10(data1[yvar])
#
#     reps = range(max(data1['rep']) + 1)
#     for rep in reps:   #とりあえず同じだし
#         data1_rep = data1.loc[data1['rep'] == rep]
#
#         slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
#         lst.append(('Monoculture', met, rep, slope_mono))
#
# file = open(f'frac_slopes.csv', 'a')  # write custom text to front
# final_pd = pd.DataFrame(lst[1:], columns=list(lst[0]))
# final_pd.to_csv('frac_slopes.csv', header=True, index=False, mode="a")
# final_pd.to_csv(f'frac_slopes.csv', index=False)



# phases = list(range(0, 20)) + [20, 30, 45]
# for phase in phases:
# #     print(phase)#
# # #     print(met)  #
# # #
# #     run(seed, "mono", 10, (0.01, 0), 10, (met, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0))
# #     run(seed, "co", 5, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1))
# #     run(seed, "mono", 10, (0.01, 0), 10, (1000, 1000, 0), (10, 0), (1, 0), 5, (3, 0), 42, "null", (1.1, 0), phase)
# #     run(seed, "co", 10, (0.01, 0.01), 10, (1, 1000, 0), (5, 5), (1, 1), 5, (3, 3), 42, "null", (1.1, 1.1), phase)
#     #     run_calc_frac(f'hcb_sim_mono_{seed}_met{met}.csv', 1000, 5)
#     run_calc_frac(f'hcb_sim_mono_685_met1000_phase0{phase}.csv', 1000, 5)
#     #     #run_calc_frac(f'hcb_sim_co_{seed}_met1.csv', 1000, 5)
#     run_calc_frac(f'hcb_sim_co_685_met1_phase0{phase}.csv', 1000, 5)
#     #
#     # compute per rep
#     data1 = pd.read_csv(f'frac_hcb_sim_mono_685_met1000_phase0{phase}.csv', header=2, na_filter=False)
#     data2 = pd.read_csv(f'frac_hcb_sim_co_685_met1_phase0{phase}.csv', header=2, na_filter=False)
#
#     yvar = 'frac'
#     data1[f'log10_{yvar}'] = np.log10(data1[yvar])
#     data2[f'log10_{yvar}'] = np.log10(data2[yvar])
#
#     reps = range(max(data1['rep']) + 1)
#     #reps2 = range(max(data2['rep']) + 1)
#     for rep in reps:   #とりあえず同じだし
#         data1_rep = data1.loc[data1['rep'] == rep]
#         data2_rep = data2.loc[data2['rep'] == rep]
#
#         slope_mono, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
#         lst.append(('Monoculture', phase, rep, slope_mono))
#
#         slope_co, intercept, r_value, p_value, std_err = linregress(data2_rep['cycle'], data2_rep[f'log10_{yvar}'])
#         lst.append(('Coculture', phase, rep, slope_co))
#
# final_pd = pd.DataFrame(lst[1:], columns=list(lst[0]))
# final_pd.to_csv(f'frac_slopes_phase0.csv', index=False)
#

# yvar = 'frac'
# data_co = pd.read_csv('frac_hcb_sim_co_561_met1.csv', header=2)
# data_co[f'log10_{yvar}'] = np.log10(data_co[yvar])
# lst = [('culture', 'met', 'rep', 'frac_slope')]
# reps = range(max(data_co['rep']) + 1)
# for rep in reps:
#     data1_rep = data_co.loc[data_co['rep'] == rep]
#
#     slope_co, intercept, r_value, p_value, std_err = linregress(data1_rep['cycle'], data1_rep[f'log10_{yvar}'])
#     lst.append(('co', 1, rep, slope_co))
#
# final_pd = pd.DataFrame(lst[1:], columns=list(lst[0]))
# final_pd.to_csv(f'frac_slopes_co.csv', index=False)


# data = pd.read_csv('phase0/frac_slopes_phase0.csv')
# data["culture"] = "Monoculture"

# for met in mets:
#     lst.append(('Avg Lin Reg Slope for Coculture at 1 Met.', met, 0, 0.11739308753272672))
# data1 = pd.DataFrame(lst[1:], columns=list(lst[0]))
#
# data2 = pd.read_csv('frac_slopes_2000.csv')
# data2["culture"] = "Monoculture"

# data = pd.concat((data2, data, data1), axis=0, ignore_index=True)


# data1 = pd.read_csv("frac_hcb_sim_co_331_met1.csv", header=2)
# data2 = pd.read_csv("frac_hcb_sim_mono_331_met1000.csv", header=2)
# data = pd.concat((data1, data2), axis=0, ignore_index=True)

data = pd.read_csv("results_29062023/5hrfrac_slope_over_met/frac_slopes_0_40_1.csv")
x = sns.relplot(data=data, x="met", y="frac_slope", hue="culture", kind="line", ci="sd", err_style="bars", alpha=0.65, palette="husl")
x.set(xlabel="Initial Methionine", ylabel="Slope of Lin. Reg. Line of Log10(Survival Fraction) by Rep", title="Impact of Initial Methionine on Monoculture Survival Fraction Evolution")
plt.subplots_adjust(top=0.95) # use a lower number to make more vertical space
plt.show()


# data1 = pd.read_csv(input(), header=2, na_filter=False)
# yvar = 'frac'
# data1[f'log10_{yvar}'] = np.log10(data1[yvar])
# slope_mono, intercept, r_value, p_value, std_err = linregress(data1['cycle'], data1[f'log10_{yvar}'])
# print(slope_mono)

# frac_hcb_sim_co_561_met1.csv
# 0.11739308753272672