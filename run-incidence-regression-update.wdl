version 1.0

task RegressionOutput {
    input {
        File barcodes_file
    }

    Int disk_size = 1 + 5 * ceil(
    size(barcodes_file, "GB")
    )

    command <<<

        python3 <<EOF
        
        from pandas import DataFrame, read_excel
        import pandas as pd
        import numpy as np
        from collections import defaultdict, Counter
        import seaborn as sns
        from matplotlib.patches import Patch
        import random
        from patsy import dmatrices
        import statsmodels.api as sm
        import statsmodels.api as sm
        import scipy
        import matplotlib.pyplot as plt
        import json
        from numpyencoder import NumpyEncoder
        import dill
        import os
        from patsy import dmatrices
        from sklearn.model_selection import KFold

        def remove_duplicate_df(df, field = 'Haplotype'):
            duplicate_rows = pd.DataFrame.duplicated(df, field, keep='first')
            return df[~duplicate_rows]
            
        def inverse_var_weight(p_array, var_array):
            var_weighted_num= np.nansum(np.divide(p_array, 
                                                    var_array, where=var_array!=0))
            var_weighted_denom = np.nansum(np.divide([1 for _ in var_array], 
                                                        var_array,where=var_array!=0))
            weighted_mean = var_weighted_num / var_weighted_denom
            
            weighted_var = 1/np.sum(1/var_array)
            weighted_std = np.sqrt(weighted_var)
            weighted_ci = (weighted_mean - 1.96 * weighted_std,
                            weighted_mean + 1.96 * weighted_std)
            
            
            return weighted_mean, weighted_var, weighted_ci

        def is_unique(s):
            a = s.to_numpy() # s.values (pandas<0.24)
            return (a[0] == a).all() #True = homozygous, #False = Heterozygous

        def return_stats(s, convert_numpy = True):
            if convert_numpy:
                a = s.to_numpy() # s.values (pandas<0.24)
            else:
                s=a
            missing_screen = a != 'X'
            poly_screen = a != 'N'
            
            if missing_screen.all() == False:
                return np.nan #missing data present in this comparison, skip
            elif poly_screen.all() == True:
                return 1. -(a[0] == a).all() #0 = homozygous, #1 = Heterozygous
                #status =  1. -(a[0] == a).all() #0 = homozygous, #1 = Heterozygous
                #if status == 1.0:
                #    random_p = np.random.random()
                #    if random_p < Ftrue:
                #        return 1. - status
                #    else:
                #        return status
                #else: 
                #    return status
                
            else: #N present in one of the comparisons, skip for now
                return np.nan

        def quantify_het(barcode_comparison, Ftrue=0, axis = 0):
            null_het= np.nansum(barcode_comparison, axis = axis) * (1-Ftrue)
            return null_het

        def quantify_total(barcode_comparison, axis = 0):
            return np.sum(~np.isnan(barcode_comparison), axis = axis)
            
        def run_fn_tests():
            test_het_calculator()
            test_unique()
            
        def test_het_calculator():
            expected_het, expected_total = (12,16)
            test_barcode1 = 6 * ['A'] + 6*['T'] + 6* ['C'] + 6 * ['G'] + 6*['X'] + 6*['N']
            test_barcode2 = 6 * ['A', 'T', 'C', 'G', 'X', 'N']
            test_df = DataFrame([test_barcode1, test_barcode2])
            b = test_df.apply(return_stats, axis=0).to_list()
            heterozygotes = np.sum(quantify_het(b), axis = 0)
            total = np.sum(quantify_total(b), axis = 0)
            assert(heterozygotes == expected_het), "Heterozygote sites should be {expected_het}, identified {x}".format(expected_het = expected_het, x = heterozygotes)
            assert(total == 16), "Total usable sites should be {expected_total}, identified {x}".format(x = total, expected_total = expected_total)

        def test_unique():
            test_df1 = DataFrame([['A'], ['B']])
            test_df2 = DataFrame(['A'], ['A'])
            assert(is_unique(DataFrame(test_df1)) == False), 'Failed unique test, misidentified [A,B] as being composed of a single unique entity'
            assert(is_unique(DataFrame(test_df2)) == True), 'Failed unique test, misidentified [A,A] as being composed of a multiple unique entities'

        def check_combinations(combo_barcodes):
            arr_2d = combo_barcodes.to_numpy()
            assert(np.array_equal(arr_2d[0], arr_2d[1]) == False), 'Duplicate barcodes were sampled'

        def calculate_RH_sample(H_mono, H_poly):
            if H_mono > H_poly:
                RH = (H_mono - H_poly) / H_mono
            else:
                RH=(H_mono - H_poly) / H_mono
            return RH

        def bootstrap(array, iterations):
            resample = np.random.choice(array, (iterations, len(array)), replace=True)
            return np.mean(resample, axis = 1)

        def RH_classifier(RH):
            RH = float(RH)
            if RH > 0.4:
                return 'cotx'
            elif RH > 0.3:
                return 'cotx_probable'
            elif RH > -0.1:
                return 'coi=2'
            elif RH > -0.2:
                return 'coi=2_probable'
            elif RH > -0.4:
                return 'coi=3'
            elif RH > -0.6:
                return 'coi=3_probable'
            elif RH > -0.8:
                return 'coi=4_probable'
            else:
                return 'coi=4'
            
        def calculate_Shannons_index(p_array):
            H = -np.sum([p_array * np.log(p_array)])
            return H

        def calculate_Evenness(H_array, k):
            E = H_array / np.log(k)
            return E

        def calculate_H12(p_array):
            p_array = sorted(p_array)
            p1 = p_array[0]
            p2 = p_array[1]
            other_p = np.asarray(p_array[2:])
            sum_rest = np.sum(other_p**2)
            H12 = (p1 + p2)**2 + sum_rest
            return H12

        def interpret_cotx_simbinbarcode(sim, cotx_event = 1):
            '''wokrks only for coi =2'''
            sim = sim[str(cotx_event)]
            converted_barcodes = []
            for binbarcode in sim:
                if binbarcode == 0:
                    barcode = 24 * [0]
                else:
                    barcode = [int(i) for i in list('{0:0b}'.format(binbarcode))]
                    while len(barcode) < 24:
                        barcode = [0] + barcode #the conversion does not preserve leading zeros
                converted_barcodes.append(barcode)
            return np.asarray(converted_barcodes)

        def wilson(p, n, z = 1.96):
            denominator = 1 + z**2/n
            centre_adjusted_probability = p + z*z / (2*n)
            adjusted_standard_deviation = np.sqrt((p*(1 - p) + z*z / (4*n)) / n)
            
            lower_bound = (centre_adjusted_probability - z*adjusted_standard_deviation) / denominator
            upper_bound = (centre_adjusted_probability + z*adjusted_standard_deviation) / denominator
            return (lower_bound, upper_bound)

        def create_k_folds(data, k=5, random_state=None): 
                kf = KFold(n_splits=k, shuffle=False, random_state=random_state)
                inclusion_idxes = [data[x[0]] for x in kf.split(data)]
                return inclusion_idxes

        def create_jackknife_lists(master_df, k = 10):
            loci_position_names = list(master_df.columns[7:31] )
            scramble_idxes,k_fold_idxes = {},{}
            for year, year_df in master_df.groupby('Year'):
                #print(year, year_df)
                scramble_idxes[year] = np.random.permutation(list(year_df.index))
                k_fold_idxes[year] = create_k_folds(scramble_idxes[year], k=k)
            
            scramble=defaultdict(list)
            for k_idx in range(k):
                for year in k_fold_idxes:
                    scramble[k_idx+1] += list(k_fold_idxes[year][k_idx])
            return scramble

        class Barcode_Stats:
            cpalette_converter = {'green': "viridis",
                                'blue':'mako',
                                'cornflowerblue': 'mako',
                                    'crimson': 'rocket',
                                    'orange': 'Oranges',
                                'purple': 'Purples_r'}
            
            
            def __init__(self, file, ISO3, sheet_name=None, adjustedN= False, k_idx = None, scramble_idxes=None): #k_idx range from 0-9
                #k_validation is a number between 1 and 10 to do k-fold leave one block out validation
                #scramble_idxes is a dictionary where that contain scrambled idxes associated with the k-th fold
                self.ISO3 = ISO3
                
                if not sheet_name:
                    sheet_name = ISO3
                    self.master_df = DataFrame(read_excel(file, sheet_name = ISO3))
                else:
                    self.master_df = DataFrame(read_excel(file, sheet_name = sheet_name))
                self.loci_position_names = list(self.master_df.columns[7:31] )
                
                if k_idx:#must reset the master_df from which everything is quantified from before doing everything eles
                    self.master_df = self.master_df.iloc[scramble_idxes[k_idx]].reset_index() #reset the master_df to be the k-fold version
                self.extract_region(ISO3)
                self.extract_mono_poly_stats()
                self.calc_stats()
                self.calculate_polyhet_timeseries()
                self.quantify_observed_poly_het(adjustedN)
                self.calculate_RH()
            
            def label_haplotypes(self):
                barcode_haplotypes_dict = {}
                unique_barcode_haplotype = 0
                barcode_haplotypes = []
                for row in self.barcodes_df[self.barcodes_df.columns[7:31]].to_numpy():
                    barcode = ('').join(row)
                    if barcode not in barcode_haplotypes_dict:
                        barcode_haplotypes_dict[barcode] = unique_barcode_haplotype
                        unique_barcode_haplotype += 1
                    barcode_haplotypes.append(barcode_haplotypes_dict[barcode])
                self.barcodes_df['Haplotype'] = barcode_haplotypes
                
            
            def extract_region(self, ISO3):
                tmp_df =self.master_df[self.master_df['ISO3'] == ISO3] #Thies samples only
                for position in self.loci_position_names:
                    tmp_df[position] = [x.strip().upper() for x in tmp_df[position]] #noticed a tab formatting in these columns
                self.barcode_pos_dict = {}
                fin = open("/usr/local/bin/barcode_pos.txt")
                for i,line in enumerate(fin.readlines()[1:]):
                    line =line.strip().split('\t')
                    self.barcode_pos_dict[self.loci_position_names[i]] = line[0] + ':' + line[1]
                
                
                tmp_df['M/P'] = [x.strip() for x in tmp_df['M/P']]
                
                
                control_samples = ['3D7', '3D7-1', '3D7-2', '3d7', 'Dd2-MOD', 'Dd2-Mod', 'Dd2/Mod','Dd2_MOD', 'DD2', 'Dd2']

                self.barcodes_df = tmp_df[(tmp_df['X'] <= 2) 
                                            & (~tmp_df['Sample_Name'].isin(control_samples))]
                self.mono_barcode_df = self.barcodes_df[(self.barcodes_df['X'] <= 2) 
                                            & (self.barcodes_df['M/P'] == 'M') 
                                            & (~self.barcodes_df['Sample_Name'].isin(control_samples))]
                self.chrono_years = sorted(Counter(self.barcodes_df['Year']).keys())
                self.label_haplotypes()


            def extract_mono_poly_stats(self, fail_threshold = 2):

                self.barcode_year_dfs = {}
                self.mono_barcode_year_dfs = {}
                self.poly_barcode_year_dfs = {}
                for year in self.chrono_years:
                    self.barcode_year_dfs[year] = self.barcodes_df[(self.barcodes_df['Year'] == year) &
                                                            (self.barcodes_df['X'] <= fail_threshold)]

                    self.mono_barcode_year_dfs[year] = self.barcodes_df[(self.barcodes_df['Year'] == year) &
                                                            (self.barcodes_df['X'] <= fail_threshold) & 
                                                            (self.barcodes_df['M/P'] == 'M')]

                    self.poly_barcode_year_dfs[year] = self.barcodes_df[(self.barcodes_df['Year'] == year) &
                                                            (self.barcodes_df['X'] <= fail_threshold) & 
                                                            (self.barcodes_df['M/P'] == 'P')]
                    self.mono_barcode_year_dfs[year].reset_index(drop=True, inplace=True)
                    self.poly_barcode_year_dfs[year].reset_index(drop=True, inplace=True)
                self.n_singles = np.asarray([len(self.mono_barcode_year_dfs[year]) for year in self.chrono_years])
                self.n_poly = np.asarray([len(self.poly_barcode_year_dfs[year]) for year in self.chrono_years])
                
                self.n_total = np.asarray([len(self.barcode_year_dfs[year]) for year in self.chrono_years])

                self.loci_allele_dict = {}
                for column in self.mono_barcode_df.columns[7:31]:
                    counts = Counter(self.mono_barcode_df[column].to_numpy())
                    counts.pop('X', 0)
                    counts.pop('N', 0)
                    major_allele = max(counts, key = counts.get)
                    counts.pop(major_allele,0)
                    try:
                        minor_allele = max(counts,key=counts.get)
                    except:
                        minor_allele = None
                    self.loci_allele_dict[column] = (major_allele, minor_allele)
                
                self.haplotype_counts = {}
                for year in self.chrono_years:
                    counts = Counter(self.barcode_year_dfs[year][self.barcode_year_dfs[year]['M/P'] == 'M']['Haplotype'])
                    self.haplotype_counts[year] = counts
                    
                
            def calc_stats(self):
                self.popgen_stats = defaultdict(list)
                n = self.n_singles + self.n_poly
                p = self.n_poly / n
                q = 1. - p
                variances = p*q /n
                wilson_interval = [wilson(prop, n_samples) for prop, n_samples in zip(p,n)]

                self.popgen_stats['p_poly_fract'] = list(p)
                self.popgen_stats['var_poly_fract'] = list(variances)
                weighted_mean, weighted_var, weighted_ci = inverse_var_weight(p, p*q/n)
                self.popgen_stats['poly_fract_inv_var'] = (weighted_mean, weighted_var, weighted_ci)
                self.popgen_stats['poly_wilson'] = wilson_interval
                x = np.asarray(self.chrono_years)
                X = sm.add_constant(x)
                y = np.array(p)
                model = sm.OLS(y,X).fit()
                self.popgen_stats['poly_fract_model'] = model

                #print('Unique Mono Fract')
                self.unique_mono_count = {}
                self.repeat_haps = defaultdict(lambda: defaultdict(lambda:0))
                for year in self.chrono_years:
                    sorted_hid= sorted(self.haplotype_counts[year], key =lambda x: self.haplotype_counts[year][x], reverse = True)
                    for hid in sorted_hid:
                        if self.haplotype_counts[year][hid] != 1:
                            self.repeat_haps[year][hid] = self.haplotype_counts[year][hid]
                        else:
                            self.repeat_haps[year]['unique'] += 1
                    self.unique_mono_count[year] = self.repeat_haps[year]['unique']
                
                self.total_mono = np.asarray([np.sum(list(self.haplotype_counts[year].values())) for year in self.chrono_years])
                p_mono_unique = np.asarray([self.unique_mono_count[year] for year in self.chrono_years])/ self.total_mono
                p_mono_clonal = (1. - p_mono_unique)
                var_mono_unique = (p_mono_unique * (1. - p_mono_unique)) / self.total_mono
                self.popgen_stats['p_mono_unique'] = list(p_mono_unique)
                self.popgen_stats['var_mono_unique'] = list(var_mono_unique)
                self.popgen_stats['wilson_mono_unique'] = [wilson(p, n) for p, n in zip(p_mono_unique, self.total_mono)]
                x = np.asarray(self.chrono_years)
                X = sm.add_constant(x)
                y = np.array(p_mono_unique)
                model = sm.OLS(y,X).fit()
                self.popgen_stats['mono_unique_model'] = model 
                weighted_mean, weighted_var, weighted_ci = inverse_var_weight(p_mono_unique, var_mono_unique)
                self.popgen_stats['mono_unique_inv_var']  = (weighted_mean, weighted_var, weighted_ci)
                
                self.popgen_stats['p_mono_clonal'] = list(p_mono_clonal) 
                self.popgen_stats['var_mono_clonal'] = list(var_mono_unique)#it's the same because it is the inverse
                self.popgen_stats['wilson_mono_clonal'] = [wilson(p, n) for p, n in zip(p_mono_clonal, self.total_mono)]
                x = np.asarray(self.chrono_years)
                X = sm.add_constant(x)
                y = np.array(p_mono_clonal)
                model = sm.OLS(y,X).fit()
                self.popgen_stats['mono_clonal_model'] = model
                weighted_mean, weighted_var, weighted_ci = inverse_var_weight(p_mono_clonal, var_mono_unique)
                self.popgen_stats['mono_clonal_inv_var']  = (weighted_mean, weighted_var, weighted_ci)

                #calculate mono diversity - all
                for year in self.chrono_years:
                    hap_ids = np.asarray(list(self.haplotype_counts[year].keys()))
                    hap_counts = np.asarray(list(self.haplotype_counts[year].values()))
                    hap_freqs = hap_counts / np.sum(hap_counts)
                    sampling_idxes = np.random.choice(hap_ids, p=hap_freqs, size = (200, 200))
                    shannon_idxes, evenness_scores, H12_scores = [],[],[]
                    for sampling_idx in sampling_idxes:
                        sampled_counts = Counter(sampling_idx)
                        sampled_freqs = np.asarray(list(sampled_counts.values())) / 200
                        
                        H12 = calculate_H12(sampled_freqs)
                        shannon_idx = calculate_Shannons_index(sampled_freqs)
                        shannon_idxes.append(shannon_idx)
                        evenness = calculate_Evenness(shannon_idx, len(list(sampled_counts.keys())))
                        evenness_scores.append(evenness)
                        H12_scores.append(H12)
                        
                    self.popgen_stats['shannon_idx_mean'].append(np.mean(shannon_idxes))
                    self.popgen_stats['evenness_mean'].append(np.mean(evenness_scores))
                    self.popgen_stats['H12_mean'].append(np.mean(H12_scores))
                    
                    self.popgen_stats['shannon_idx_var'].append(np.var(shannon_idxes))
                    self.popgen_stats['evenness_var'].append(np.var(evenness_scores))
                    self.popgen_stats['H12_var'].append(np.var(H12_scores))
                    
                self.popgen_stats['shannon_idx_mean'] =  np.array(self.popgen_stats['shannon_idx_mean'])
                model = sm.OLS(self.popgen_stats['shannon_idx_mean'],X).fit()
                self.popgen_stats['shannon_idx_model'] = model
                
                self.popgen_stats['evenness_mean'] =  np.array(self.popgen_stats['evenness_mean'])
                model = sm.OLS(self.popgen_stats['evenness_mean'],X).fit()
                self.popgen_stats['evenness_model'] = model
                
                self.popgen_stats['H12_mean'] =  np.array(self.popgen_stats['H12_mean'])
                model = sm.OLS(self.popgen_stats['H12_mean'],X).fit()
                self.popgen_stats['H12_model'] = model
                
                if 'mccoil_median' in self.barcodes_df.columns:
                    for year, df in self.barcodes_df.groupby('Year'):
                        self.popgen_stats['mccoil_coi'].append(np.mean(df['mccoil_median']))
                        self.popgen_stats['mccoil_coi_std'].append(np.std(bootstrap(df['mccoil_median'], 1000)))
                    for year, df in self.barcodes_df.groupby('Year'):
                        self.popgen_stats['mccoil_coi_poly'].append(np.mean([x for x in df['mccoil_median'] if x >=2]))
                        self.popgen_stats['mccoil_coi_poly_std'].append(np.std([x for x in df['mccoil_median'] if x >=2]))
                
            def calculate_polyhet_timeseries(self):
                self.poly_het_dict = defaultdict(dict)
                for year in self.chrono_years:
                    for position in self.loci_position_names:
                        counts = Counter(self.poly_barcode_year_dfs[year][position].to_list())
                        n_missing = counts.pop('X', 0)
                        total = np.sum(list(counts.values()))
                        n_het = counts.pop('N',0)
                        p_het = n_het  / total
                        #print(year, position, p_het, n_missing, n_missing/total, total)
                        self.poly_het_dict[year][position] = (p_het, total)

                self.poly_het_timeseries = defaultdict(list)
                for loci in self.loci_position_names:
                    for year in self.chrono_years:
                        self.poly_het_timeseries[loci].append(self.poly_het_dict[year][loci][0])
            
            def quantify_observed_poly_het(self, adjustedN = False):
                self.poly_barcode_het_dist = defaultdict(list)
                self.poly_barcode_het_avg = {}
                self.poly_samples = {}
                if not adjustedN:
                    #counting N from barcode
                    for year in self.chrono_years:
                        for i,row in enumerate(self.poly_barcode_year_dfs[year][self.poly_barcode_year_dfs[year].columns[7:31]].to_numpy()):
                            barcode = ('').join(row)
                            barcode_counts = Counter(barcode)
                            barcode_counts.pop('X',0)
                            total = np.sum(list(barcode_counts.values()))
                            het = barcode_counts.pop('N')
                            H_poly_barcode = het/total
                            self.poly_barcode_het_dist[year].append(H_poly_barcode)
                        self.poly_samples[year] = list(self.poly_barcode_year_dfs[year]['Sample_Name'])
                        self.poly_barcode_het_avg[year] = np.mean(self.poly_barcode_het_dist[year])
                else:
                    for year in self.chrono_years:
                        total = 24. - np.asarray(self.poly_barcode_year_dfs[year]['X'])
                        het = np.asarray(self.poly_barcode_year_dfs[year]['Adjusted_Het']) / total
                        self.poly_barcode_het_dist[year] = het
                        self.poly_samples[year] = list(self.poly_barcode_year_dfs[year]['Sample_Name'])
                        self.poly_barcode_het_avg[year] = np.mean(het)
                
            def sample_poly_barcodes(self,year, coi=2, samples = 100, Ftrue = 0):
                combinations = np.random.randint(self.mono_barcode_year_dfs[year].index[0], 
                                                self.mono_barcode_year_dfs[year].index[-1], 
                                                size = (samples,coi))
                sampled_poly_barcodes = []
                for i, combo in enumerate(combinations):
                    sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[combo, ['Haplotype']]
                    flag = len(np.unique(list(sampled_haplotypes['Haplotype']))) !=coi
                    while flag:
                        combo = np.random.randint(self.mono_barcode_year_dfs[year].index[0], 
                                                self.mono_barcode_year_dfs[year].index[-1], coi)
                        combinations[i] = combo
                        sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[combo, ['Haplotype']]
                        flag = len(np.unique(list(sampled_haplotypes['Haplotype']))) !=coi
                    combo_barcodes = self.mono_barcode_year_dfs[year].loc[combo, self.loci_position_names]
                    check_combinations(combo_barcodes)
                    #sampled_poly_barcodes.append(combo_barcodes.apply(lambda x: return_stats(x, Ftrue), axis = 0).to_list())
                    sampled_poly_barcodes.append(combo_barcodes.apply(lambda x: return_stats(x), axis = 0).to_list())

                sampled_poly_barcodes = np.asarray(sampled_poly_barcodes)

                return sampled_poly_barcodes

            def simulator(self,n_poly, n_iterations = 1000, Ftrue = 0, axis = 0, coi = 2):
                run_fn_tests()
                poly_simulations = defaultdict(list)
                for year, n in zip(self.chrono_years, n_poly):
                    for n_rep in range(n_iterations):
                        attempt_count = 1
                        #if n_rep % 100 ==0:
                        #    print(year,n, n_rep)
                        b = self.sample_poly_barcodes(year,coi=coi, samples=n)
                        heterozygotes = quantify_het(b, Ftrue, axis = axis)
                        total = quantify_total(b, axis = axis)  
                        p_het = heterozygotes / total
                        while (np.isfinite(p_het).all() == False): #if zero is found in the total
                            assert(attempt_count <= 3), 'maximum number of attempts reached'
                            #print('attempting resample {x}'.format(x=attempt_count))
                            b = self.sample_poly_barcodes(year, samples=n)
                            heterozygotes = quantify_het(b, Ftrue, axis = axis)
                            total = quantify_total(b, axis = axis)
                            p_het = heterozygotes/ total
                            attempt_count += 1
                        poly_simulations[year].append(p_het)
                    poly_simulations[year] = np.asarray(poly_simulations[year])
                return poly_simulations
            
            def sample_cotx_barcodes_from_mono(self, year, coi=2, samples = 100, initial_coi = 2, cotx_event = 1):
                #cotx_sim stats
                cotx_sim_filtered_data = [x for x in self.cotx_simulation_data[initial_coi]['cotx_event{x}'.format(x=cotx_event)] if len(x) >=coi]
                random_sim_idxes = np.asarray(random.sample(range(len(cotx_sim_filtered_data)), samples))
                sampled_cotx_barcodes = np.asarray(cotx_sim_filtered_data, dtype ='object')[random_sim_idxes]

                #superinfection layer
                combinations = np.random.randint(self.mono_barcode_year_dfs[year].index[0], 
                                                self.mono_barcode_year_dfs[year].index[-1], 
                                                size = (samples,initial_coi))
                sampled_poly_barcodes = []
                for i, combo in enumerate(combinations):
                    sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[combo, ['Haplotype']]
                    flag = len(np.unique(list(sampled_haplotypes['Haplotype']))) !=initial_coi
                    while flag:
                        combo = np.random.randint(self.mono_barcode_year_dfs[year].index[0], 
                                                    self.mono_barcode_year_dfs[year].index[-1], initial_coi)
                        combinations[i] = combo
                        sampled_haplotypes = self.mono_barcode_year_dfs[year].loc[combo, ['Haplotype']]
                        flag = len(np.unique(list(sampled_haplotypes['Haplotype']))) !=initial_coi
                    combo_barcodes = self.mono_barcode_year_dfs[year].loc[combo, self.loci_position_names]
                    #print(combo, combo_barcodes)
                    check_combinations(combo_barcodes)

                    #print(sampled_cotx_barcodes[i])
                    relatedness_maps = []
                    for n in range(coi):
                        relatedness_maps.append(sampled_cotx_barcodes[i][n])
                        #relatedness_map1 = sampled_cotx_barcodes[i][0]
                        #relatedness_map2 = sampled_cotx_barcodes[i][1]

                    combo_barcodes = combo_barcodes.to_numpy()
                    #print(relatedness_maps)
                    #print(combo_barcodes)
                    cotx_strains = []
                    for relatedness_map in relatedness_maps:
                        tmp = [combo_barcodes[strain_choice][position] for position,strain_choice in enumerate(relatedness_map)]
                        cotx_strains.append(tmp)
                    #cotx_strain1 = [combo_barcodes[strain_choice][position] for position,strain_choice in enumerate(relatedness_map1)]
                    #cotx_strain2 = [combo_barcodes[strain_choice][position] for position,strain_choice in enumerate(relatedness_map2)]
                    #print('Cotx Strains')
                    #print(cotx_strains)
                    df = DataFrame(cotx_strains)#[cotx_strain1, cotx_strain2])
                    sampled_poly_barcodes.append(df.apply(lambda x: return_stats(x), axis = 0).to_list())



                sampled_poly_barcodes = np.asarray(sampled_poly_barcodes)


                return sampled_poly_barcodes

            def simulator_cotx(self, n_poly, n_iterations = 1000, axis = 0, coi = 2, initial_coi = 2, cotx_event=1):
                run_fn_tests()
                poly_simulations = defaultdict(list)
                for year, n in zip(self.chrono_years, n_poly):
                    for n_rep in range(n_iterations):
                        attempt_count = 1
                        #if n_rep % 100 ==0:
                        #    print(year,n, n_rep)
                        b = self.sample_cotx_barcodes_from_mono(year,coi=coi, samples=n, initial_coi = initial_coi,cotx_event=cotx_event)
                        heterozygotes = quantify_het(b, Ftrue=0, axis = axis)
                        total = quantify_total(b, axis = axis)  
                        p_het = heterozygotes / total
                        while (np.isfinite(p_het).all() == False): #if zero is found in the total
                            assert(attempt_count <= 3), 'maximum number of attempts reached'
                            #print('attempting resample {x}'.format(x=attempt_count))
                            b = self.sample_poly_barcodes(year, samples=n)
                            heterozygotes = quantify_het(b, Ftrue, axis = axis)
                            total = quantify_total(b, axis = axis)
                            p_het = heterozygotes/ total
                            attempt_count += 1
                        poly_simulations[year].append(p_het)
                    poly_simulations[year] = np.asarray(poly_simulations[year])
                return poly_simulations

            def calculate_RH(self, n_poly_per_year=200, n_iter=200): #can set to 20 for fast
                self.RH_barcode_dict  = {}
                self.observed_RH = defaultdict(list)
                
                #print('Simulating mono barcode sampling for RH')
                H_mono_barcodes = self.simulator([n_poly_per_year for x in self.n_poly], n_iter, 0, axis = 1)
                self.calculate_RH_barcode_distribution(H_mono_barcodes)
                self.calculate_RHyear_distribution(H_mono_barcodes)
                self.H_mono_barcodes=H_mono_barcodes
                
                
                #actual samples
                minimum_sample = np.min(self.n_poly)
                for year in self.chrono_years:
                    for i, H_poly_barcode in enumerate(self.poly_barcode_het_dist[year]):
                        RH_sample_dist = [calculate_RH_sample(np.mean(sim_trace), H_poly_barcode) for sim_trace in self.H_mono_barcodes[year]]
                        self.observed_RH[year].append(RH_sample_dist)
                        
                X = sm.add_constant(self.chrono_years)
                y = np.array([np.mean(self.RH_barcode_dict[year]) for year in self.chrono_years])
                model1 = sm.OLS(y,X).fit()
                self.RH_barcode_dict['model'] = model1
                
                data,Rh_sample_averages,classification_averages = [], defaultdict(list), defaultdict(lambda: defaultdict(list))
                for year in self.chrono_years:
                    for sample_name, RH_sample_dist in zip(self.poly_samples[year], self.observed_RH[year]):
                        data.append([sample_name, year, np.mean(RH_sample_dist)])
                self.RH_df = DataFrame(data)
                self.RH_df.columns=['Sample', 'Year', 'RH']
                self.RH_df['classification'] = self.RH_df['RH'].apply(RH_classifier)
                self.poly_df = pd.merge(self.master_df, self.RH_df, left_on = 'Sample_Name', right_on = 'Sample')
                
                
                for year, df in self.RH_df.groupby('Year'):
                    total_sample =np.asarray(df['RH'])
                    cotx_total_samples = np.asarray(df['classification'])
                    
                    sampling_idxes = np.random.randint(0,len(total_sample), size = (100,200))
                    for sampling_idx in sampling_idxes:
                        RH_sample = total_sample[sampling_idx]
                        Rh_sample_averages[year].append(np.mean(RH_sample))
                        
                        cotx_samples = cotx_total_samples[sampling_idx]
                        cotx_counts = Counter(cotx_samples)
                        p_cotx =(cotx_counts['cotx'] + cotx_counts['cotx_probable'])  / len(sampling_idx)
                        classification_averages['cotx'][year].append(p_cotx)
                        
                        for classification in ['coi=2', 'coi=3', 'coi=4']:
                            p =(cotx_counts[classification] + cotx_counts[classification + '_probable'])  / len(sampling_idx)
                            if classification != 'cotx':
                                classification_averages[classification][year].append(p)
                            
                averages =  np.asarray([np.mean(Rh_sample_averages[year]) for year in self.chrono_years])
                variances = np.asarray([np.var(Rh_sample_averages[year]) for year in self.chrono_years])
                weighted_mean, weighted_var, weighted_ci = inverse_var_weight(averages, variances)
                
                
                self.RH_yearly_averages =averages
                self.RH_yearly_variances = variances
                self.RH_average = weighted_mean
                self.RH_weighted_var = weighted_var
                self.RH_ci = weighted_ci
                
                for classification in classification_averages:
                    averages =  np.asarray([np.mean(classification_averages[classification][year]) for year in self.chrono_years])
                    variances = np.asarray([np.var(classification_averages[classification][year]) for year in self.chrono_years])
                    weighted_mean, weighted_var, weighted_ci = inverse_var_weight(averages, variances)
                    self.popgen_stats[classification + '_average'] = averages
                    self.popgen_stats[classification + '_var'] = variances
                    self.popgen_stats[classification + '_inv_var'] = (weighted_mean, weighted_var, weighted_ci)
                
            
            
            def calculate_RH_barcode_distribution(self, poly_simulations):
                def distance_RH_trace(RH, H_mono_trace, H_poly_trace):
                    H_mono_trace = np.asarray(H_mono_trace)
                    H_poly_trace = np.asarray(H_poly_trace)
                    distance = np.sum((H_mono_trace - (H_mono_trace * RH) - H_poly_trace)**2)
                    return distance
            
                H_poly_trace = list(self.poly_barcode_het_avg.values())
                n_reps = len(poly_simulations[list(poly_simulations.keys())[0]])
                barcode_het_timetraces = []
                for irep in range(n_reps):
                    timetrace= []
                    for year in self.chrono_years:
                        timetrace.append(np.mean(poly_simulations[year][irep]))
                    barcode_het_timetraces.append(timetrace)

                RH_barcode_distribution= []
                for itrace in range(n_reps):
                    output = scipy.optimize.minimize(lambda RH: distance_RH_trace(RH, barcode_het_timetraces[itrace],H_poly_trace),
                                                    x0 = [0.5], bounds = [(-1,1)])
                    RH_barcode_distribution.append(output)
                self.RH_barcode_dict['total'] = [output.x[0] for output in RH_barcode_distribution]

            def calculate_RHyear_distribution(self, poly_simulations):
                def distance_RHyear_trace(RH, H_mono_traces_dict, itrace, year, H_poly_dict = self.poly_barcode_het_dist):
                    distance = 0
                    H_mono_trace = np.mean(H_mono_traces_dict[year][itrace])
                    H_poly_trace = np.mean(H_poly_dict[year])
                    distance_array = (H_mono_trace - (H_mono_trace* RH) - H_poly_trace) **2
                    return np.sum(distance_array)
            
                for year in self.chrono_years:
                    n_reps = len(poly_simulations[list(poly_simulations.keys())[0]])
                    RH_year_distribution = []
                    for itrace in range(n_reps):
                        output = scipy.optimize.minimize(lambda RH: distance_RHyear_trace(RH, poly_simulations, itrace, year),
                                                        x0 = [0.5], bounds = [(0,1)])
                        RH_year_distribution.append(output)

                    self.RH_barcode_dict[year] =[output.x[0] for output in RH_year_distribution]
            
            def simulate_coi_cotx_sweep(self,n_poly = 200, n_iter = 200, oocyst_alpha = 2.5):
                self.model_expectations = defaultdict(lambda: defaultdict(list))
                self.cotx_simulation_data = load_cotx_simulations(oocyst_alpha)
                
                H_mono_barcode_coi=defaultdict(dict)
                for coi in [2, 3,4,5]:
                    print('Simulating COI={coi} expectation'.format(coi=coi))
                    H_mono_barcode_coi[coi] = self.simulator([n_poly for x in self.n_poly], n_iter, 0,  
                                                    axis = 1, coi = coi)

                for year in self.chrono_years:
                    for coi in [2,3,4,5]:
                        for sim_iteration in H_mono_barcode_coi[coi][year]:
                            mean_het = np.mean(self.H_mono_barcodes[year])
                            for sim in sim_iteration:
                                self.model_expectations['coi={x}, oocyst_alpha={alpha}'.format(x=coi,alpha = oocyst_alpha)][year].append(calculate_RH_sample(mean_het, sim))
                
                for initial_coi in [2,3,4,5]:
                    for cotx_round in [1,2,3]:
                        print('Initial COI = {initial_coi}, Simulating cotx round {x}'.format(initial_coi = initial_coi,
                                                                                            x=cotx_round))
                        H_mono_barcode_cotx = self.simulator_cotx([n_poly for x in self.n_poly], n_iter, axis = 1, 
                                                            coi=2, initial_coi=initial_coi,cotx_event =cotx_round)
                        for year in self.chrono_years:
                            for sim_iteration in H_mono_barcode_cotx[year]:
                                mean_het = np.mean(self.H_mono_barcodes[year])
                                for sim in sim_iteration:
                                    self.model_expectations['cotx_{initial_coi}_{alpha}_{cotx_round}'.format(initial_coi=initial_coi,
                                                                                                alpha = oocyst_alpha, cotx_round = cotx_round)][year].append(calculate_RH_sample(mean_het, sim))

            
            def plot_sample_distribution(self, color, ax = None, x_annotate=-0.2, y_annotate =0.1, legend = True, title = None):
                if not ax:
                    fig, ax = plt.subplots()
                ax.bar([year for year in self.chrono_years], 
                        self.n_singles, color = 'grey',alpha = 0.3) 
                        #yerr = np.sqrt(n_total * n_singles/n_total * n_poly/n_total))

                bar = ax.bar([year for year in self.chrono_years], 
                        self.n_poly, color = color, bottom = self.n_singles)#,  yerr = np.sqrt(n_total * n_singles/n_total * n_poly/n_total))
                
                for x,y,z in zip(self.chrono_years, self.n_singles, self.n_singles):
                    x = round(x,2)
                    y = round(y,2)
                    ax.annotate(str(z),(x+x_annotate, y+y_annotate), 
                                fontsize = 15, color = 'black',fontweight="bold")
                    
                for x,y,z in zip(self.chrono_years, self.n_singles + self.n_poly, self.n_poly):
                    x = round(x,2)
                    y = round(y,2)
                    ax.annotate(str(z),(x+x_annotate, y+y_annotate), 
                                fontsize = 12, color = color,fontweight="bold")
                    
                    
                legend_elements = [Patch(facecolor='grey', edgecolor='black',label='Monogenomic'),
                                Patch(facecolor=color, edgecolor='black',label='Polygenomic')]
                if legend:
                    ax.legend(handles = legend_elements)
                if not title:
                    ax.set_title('Sample Distribution', fontsize = 20)
                else:
                    ax.set_title(title, fontsize = 20, loc = 'left')
                ax.set_xticks(self.chrono_years)
                ax.set_xticklabels(self.chrono_years)
                ax.tick_params(labelsize = 15)
                return ax
            
            def plot_mono_poly_fract(self, color, ax = None, 
                                    x_annotate = -0.02, y_annotate = 0.02, annotate_color = 'black',title = None):
                if not ax:
                    fig, ax = plt.subplots()
                p_mono = self.n_singles / self.n_total
                bar=ax.bar([year for year in self.chrono_years], p_mono, color = 'grey', alpha = 0.3)
                for b, z in zip(bar, [round(_,2) for _ in p_mono]):
                    x,y =  b._x0, b._height
                    #print(x,y)
                    if str(x) == 'nan':
                        x = 0
                    if str(y) == 'nan':
                        y = 0
                    x = round(x,2)
                    y = round(y,2)
                    ax.annotate(str(z),(x+x_annotate, y+y_annotate), fontsize = 12, 
                                color = annotate_color, fontweight="bold")
                
                ax.bar([year for year in self.chrono_years],
                    self.n_poly / self.n_total, color = color, bottom = self.n_singles / self.n_total)
                if not title:
                    ax.set_title('Mono vs Poly Fraction', fontsize = 20)
                else:
                    ax.set_title(title, fontsize = 20, loc = 'left')
                ax.tick_params(labelsize = 15)
                ax.set_xlim(self.chrono_years[0]-1,self.chrono_years[-1]+1)
                ax.set_xticks(self.chrono_years)
                ax.set_xticklabels(self.chrono_years)

            
            def plot_mono_hap_sharing(self, color, ax = None, annotate_color = 'black', x_annotate = 0.2, y_annotate = 0.05, title = None):
                if not ax:
                    fig, ax = plt.subplots()

                for year in self.chrono_years:
                    bottom, unique = 0,0
                    sorted_hid= sorted(self.haplotype_counts[year], 
                                    key =lambda x: self.haplotype_counts[year][x], 
                                    reverse = True)
                    cpalette = sns.color_palette(Barcode_Stats.cpalette_converter[color],  len(self.repeat_haps[year]))
                    total = np.sum(list(self.repeat_haps[year].values()))
                    for i,hid in enumerate(self.repeat_haps[year]):
                        if hid != 'unique':
                            height = self.repeat_haps[year][hid]/total
                            ax.bar([year], height, bottom = bottom, color = cpalette[i], edgecolor = 'black')
                            bottom += height
                    #shared_fracts.append(round(bottom,2))
                    bar = ax.bar([year], self.repeat_haps[year]['unique']/total, bottom = bottom, color = 'grey', alpha = 0.3)
                    
                shared_fracts = [round(_,2) for _ in 1. - np.asarray(self.popgen_stats['p_mono_unique'])]
                for x,y,z in zip(self.chrono_years, shared_fracts, shared_fracts):
                    x = round(x,2)
                    y = round(y,2)
                    ax.annotate(str(z),(x-x_annotate, y+y_annotate), 
                                fontsize = 12, color = annotate_color,fontweight="bold")
                    
                if not title:
                    ax.set_title('Mono Clonality', fontsize = 20)
                else:
                    ax.set_title(title, fontsize = 20, loc = 'left')
                ax.tick_params(labelsize = 15)
                ax.set_xlim(self.chrono_years[0]-1,self.chrono_years[-1]+1)
                ax.set_xticks(self.chrono_years)
                ax.set_xticklabels(self.chrono_years)
                    
            def plot_persistent_clones(self, color, ax = None, x_annotate = [-0.3, 0.2], y_annotate = 0.1, title =None):
                if not ax:
                    fig, ax = plt.subplots()
                cpalette = sns.color_palette(Barcode_Stats.cpalette_converter[color],  10)
                hap_stats = defaultdict(dict)
                total_clusters = 1
                for year in self.chrono_years:
                    sorted_hid= sorted(self.haplotype_counts[year], key =lambda x: self.haplotype_counts[year][x], reverse = True)
                    for hid in sorted_hid:
                        if self.haplotype_counts[year][hid] != 1:
                            hap_stats[hid][year] = self.haplotype_counts[year][hid]
                            total_clusters += 1
                cpalette = sns.color_palette(Barcode_Stats.cpalette_converter[color],  total_clusters)
                
                count = 1
                flipper = 0
                for haplotype in hap_stats:
                        x_array, y_array, s_array = [], [], []
                        for year in hap_stats[haplotype]:
                            x_array.append(year)
                            y_array.append(count)
                            s_array.append(hap_stats[haplotype][year] * 100)
                        ax.scatter(x_array,y_array, s_array, color = cpalette[count-1], edgecolor = 'black')
                        ax.plot(x_array,y_array, color = cpalette[count-1])
                        for i, txt in enumerate(s_array):
                            ax.annotate(int(txt/100.), 
                                        (x_array[i]+x_annotate[flipper], y_array[i]-y_annotate), 
                                        fontsize = 12, color='black', fontweight="bold")
                        if flipper == 0:
                            flipper = 1
                        else:
                            flipper =0 
                        count += 1
                if not title:
                    ax.set_title('Mono Clonality', fontsize = 20)
                else:
                    ax.set_title(title, fontsize = 20, loc = 'left')
                ax.tick_params(labelsize = 15)
                ax.tick_params(left=False,
                    labelleft=False)
                
                ax.set_xlim(self.chrono_years[0]-1,self.chrono_years[-1]+1)
                ax.set_xticks(self.chrono_years)
                ax.set_xticklabels(self.chrono_years)

            def plot_longitudinal(self, field,color='orange', ax=None, inverse_var = False, title = None):
                if not ax:
                    fig, ax = plt.subplots()
                fields = {'mono_unique': ('p_mono_unique', 'var_mono_unique', 'mono_unique_model', 'mono_unique_inv_var'),
                        'poly_fract': ('p_poly_fract', 'var_poly_fract', 'poly_fract_model','poly_fract_inv_var'),
                        'cotx': ('cotx_average', 'cotx_var', 'cotx_inv_var')}
                        #'evenness': ('evenness_mean', 'evenness_var', 'evenness_model'),
                        #'H12': ('H12_mean', 'H12_var', 'H12_model')}
            
                p = self.popgen_stats[fields[field][0]]
                variances = self.popgen_stats[fields[field][1]]
                
                if not inverse_var:
                    model = self.popgen_stats[fields[field][2]]
                    x = np.asarray(self.chrono_years)
                    X = sm.add_constant(x)
                    ax.scatter(self.chrono_years, p, color = color)
                    ax.plot(x, model.predict(X), color = color)
                    ax.fill_between(self.chrono_years, p + 2.5*np.sqrt(variances), 
                                    p-2.5*np.sqrt(variances), 
                                    alpha = 0.3, color = color)
                    ax.set_ylim(0,1)
                
                    ax.set_title(field, fontsize = 20)
                    #ax.set_xlabel('Year', fontsize = 15)
                    ax.set_ylabel('Proportion', fontsize = 15)
                    ax.tick_params(labelsize = 15)
                    ax.set_xlim(self.chrono_years[0]-1,self.chrono_years[-1]+1)
                    ax.set_xticks(self.chrono_years)
                    ax.set_xticklabels(self.chrono_years)
                
                else:
                    shifted_x = [x+1 for x in range(len(self.chrono_years))]
                    ax.errorbar(shifted_x, p, color = color, 
                    yerr = 1.96*np.sqrt(variances), 
                    markersize=20,
                    markeredgecolor=color,
                    markerfacecolor=color,
                    fmt = '.', ecolor = color, capsize = 10)
                    
                    weighted_mean, weighted_var, weighted_coi = self.popgen_stats[fields[field][3]]
                    x = [shifted_x[0], shifted_x[-1]]
                    ax.plot(x,[weighted_mean for _ in x], color = color, linewidth = 3)
                    ax.fill_between(x, [weighted_coi[0] for _ in x],
                                    [weighted_coi[1] for _ in x], 
                                    color = color, linewidth = 3, alpha = 0.2)
                    ax.set_ylim(0,1)
                    if not title:
                        ax.set_title(field, fontsize = 20)
                    else:
                        ax.set_title(title, fontsize = 20, loc = 'left')
                    #ax.set_xlabel('Year', fontsize = 15)
                    ax.set_ylabel('Proportion', fontsize = 15)
                    ax.tick_params(labelsize = 15)
                    ax.set_xlim(shifted_x[0]-0.5,shifted_x[-1]+0.5)
                    ax.set_xticks(shifted_x)
                    ax.set_xticklabels(shifted_x)
                
                
            def plot_RH_average_confidence(self, color = 'orange', ax = None):
                if not ax:
                    fig, ax = plt.subplots()
                RH = np.mean(self.RH_barcode_dict['total'])
                RH_ci = (np.percentile(self.RH_barcode_dict['total'], 2.5),
                        np.percentile(self.RH_barcode_dict['total'], 97.5))
                ax.hist(self.RH_barcode_dict['total'], color = 'orange')
                
                ax.set_xlabel(r'$R_{H}$', fontsize = 15)
                ax.set_ylabel('Freq', fontsize = 15)

            
            def plot_RHsample_longitudinal_average(self,color='orange', ax = None):
                if not ax:
                    fig, ax = plt.subplots()
                y = np.array([np.mean(self.RH_barcode_dict[year]) for year in self.chrono_years])
                X = sm.add_constant(self.chrono_years)
                
                ax.scatter(self.chrono_years, y, color = color)
                ax.plot(self.chrono_years, self.RH_barcode_dict['model'].predict(X), color = color)
                ax.boxplot([self.RH_barcode_dict[year] for year in self.chrono_years], 
                        positions = self.chrono_years, showfliers=False,
                        notch=True, patch_artist=True,
                            boxprops=dict(facecolor=color, color=color),
                            capprops=dict(color=color),
                            whiskerprops=dict(color=color),
                            flierprops=dict(color=color, markeredgecolor=color),
                            medianprops=dict(color=color),
                            )

                ax.set_ylim(0,0.5)
                ax.tick_params(axis='both', labelsize = 15)
                ax.set_ylabel(r'$R_{H}$', fontsize = 20)
                ax.set_xlabel('Year', fontsize = 20)
                
            def plot_RHsample_longitudinal(self,color = 'orange', ax=None):
                if not ax:
                    fig, ax = plt.subplots(figsize=(12,5))
                b = sns.swarmplot(x = 'Year', y= 'RH', 
                                data = self.RH_df, color = color, ax = ax)
                ax.plot([0-0.5,
                        len(self.chrono_years) + 0.5],
                        [self.RH_average, self.RH_average], 
                        color = 'black', linewidth = 3)
                ax.tick_params(axis='both', labelsize=15)
                
                sns.boxplot(showmeans=True,
                            #meanline=True,
                            meanprops={"marker":"s",
                                    "markerfacecolor":"white", 
                                    "markeredgecolor":color,"markersize":'12'},
                            medianprops={'visible': False},
                            whiskerprops={'visible': False},
                            zorder=10,
                            x="Year",
                            y='RH',
                            data=self.RH_df,
                            showfliers=False,
                            showbox=False,
                            showcaps=False,
                            ax=ax)
                ax.set_xlabel('Year', fontsize = 20)
                ax.set_ylabel(r'$R_{H}$', fontsize = 20)
                legend_elements = [Patch(facecolor='black', edgecolor='black',
                                        label=r'$R_{H}=$' + str(round(self.RH_average,2)) + 
                                        ' ' + '({ci1},{ci2})'.format(ci1 = str(round(self.RH_ci[0],2)),
                                                                    ci2 = str(round(self.RH_ci[1],2))))]
                
                ax.legend(handles = legend_elements, fontsize = 15)
                ax.set_ylim(-1.1, 1.0)
                
            def plot_cotx_sweep(self, ax = None):
                if not ax:
                    fig, ax = plt.subplots(figsize=(12,5))
                    
                simulation_boxplot_results = []
                for year in self.chrono_years:
                    for key in self.model_expectations:
                        for RH in self.model_expectations[key][year]:
                            simulation_boxplot_results.append([year,RH, key])
                df = DataFrame(simulation_boxplot_results)
                df.columns=['Year', 'RH', 'Condition']
                df_melt = df.melt(id_vars = ['Year', 'Condition'], value_vars = 'RH')
                cotx_colors = sns.color_palette('rocket', 4)
                superinfection_colors = sns.color_palette('mako_r', 3)
                order = ['cotx_5_2_3','cotx_4_2_3','cotx_3_2_3', 'cotx_2_2_3', 
                                        'cotx_5_2_2','cotx_4_2_2','cotx_3_2_2','cotx_2_2_2',
                                        'cotx_5_2_1','cotx_4_2_1','cotx_3_2_1', 'cotx_2_2_1',
                                        'coi=2', 'coi=3', 'coi=4']
                colors = 3*list(cotx_colors.as_hex()) + list(superinfection_colors.as_hex())
                custom_pal = {}
                for x,y in zip(order, colors):
                    custom_pal[x] = y

                ax.fill_between([-0.5,3.5],[1.0,1.0],[-1.5,-1.5], color = 'grey', alpha = 0.1)
                ax.fill_between([8-0.5,11.5],[1.0,1.0],[-1.5,-1.5], color = 'grey', alpha = 0.1)

                b = sns.boxplot(data = df_melt,
                                x = 'Condition',
                                y = 'value',showfliers = False, 
                                ax = ax, palette = custom_pal,
                                showmeans=True, 
                                meanprops={"marker":"o",
                                    "markerfacecolor":"white", 
                                    "markeredgecolor":"black",
                                    "markersize":"10"},
                                order= order)

                ax.tick_params(axis='both', labelsize=15)
                ax.set_ylabel(r'$R_{H}$', fontsize = 20)
                ax.set_xlabel('Condition', fontsize = 20)


                line3 = 4*['3'] + 4*['2'] + 4*['1'] + ['COI=2', 'COI=3', 'COI=4']
                ax.set_xticklabels(line3, rotation=45)

                legend_elements = [Patch(facecolor=cotx_colors[0], edgecolor='black',label=r'$COI_{i,cotx}=5$'),
                                Patch(facecolor=cotx_colors[1], edgecolor='black',label=r'$COI_{i,cotx}=4$'),
                                Patch(facecolor=cotx_colors[2], edgecolor='black',label=r'$COI_{i,cotx}=3$'),
                                Patch(facecolor=cotx_colors[3], edgecolor='black',label=r'$COI_{i,cotx}=2$'),
                                Patch(facecolor=superinfection_colors[0], edgecolor='black',label=r'$COI_{i,super}=2$'),
                                Patch(facecolor=superinfection_colors[1], edgecolor='black',label=r'$COI_{i,super}=3$'),
                                Patch(facecolor=superinfection_colors[2], edgecolor='black',label=r'$COI_{i,super}=4$')]

                ax.legend(handles = legend_elements, fontsize = 12)

                ax.plot([0,14], [0,0], color = 'black', linestyle= '--')
                ax.plot([0,14], [0.3,0.3], linestyle= '--', color = 'crimson')
                ax.annotate('Cotx Detection\n     Threshold', [12.5,0.35], color = 'crimson')
                

            def plot_RH_classification(self, color = 'orange', ax = None):
                if not ax:
                    fig, ax = plt.subplots(figsize=(12,5))
                classification_counts = Counter(self.RH_df['classification'])
                total = np.sum(list(classification_counts.values()))
                x_array = np.asarray([0,1,2,3])
                #colors = sns.color_palette(Barcode_Stats.cpalette_converter[color],3)

                cat1 =[classification_counts[key]/total for key in ['cotx', 'coi=2', 'coi=3', 'coi=4']]
                cat2 =[classification_counts[key]/total for key in ['cotx_probable', 'coi=2_probable', 'coi=3_probable', 'coi=4_probable']]
                ax.bar(x_array, cat1, color=color)
                ax.bar(x_array, cat2, bottom=cat1, color = color)

                proportions = np.asarray(cat1) + np.asarray(cat2)
                
                stdev_array = []
                for p in proportions:
                    wilson_low, wilson_high = wilson(p, total)
                    rel_low_boundary = p - wilson_low
                    rel_high_boundary = wilson_high - p
                    stdev_array.append((rel_low_boundary, rel_high_boundary))
                stdev_array = np.asarray(stdev_array).T

                ax.errorbar(x_array, proportions, yerr = stdev_array, fmt='.', capsize = 5, color = 'black')


                ax.set_xticks(x_array)
                ax.set_xticklabels(['Cotransmission', 'COI=2', 'COI=3', 'COI=4'], fontsize = 15, rotation = 45)
                ax.set_ylabel('Proportion', fontsize = 15)
                ax.tick_params(labelsize = 15)
                ax.set_ylim(0,1)
                #ax.set_title('2020', fontsize = 20)

            def generate_summary_report(self, output_file = None):
                data_report = {}
                data_report['n_singles'] = self.n_singles
                data_report['n_poly'] = self.n_poly
                data_report['n_total'] = self.n_singles + self.n_poly

                data_report['poly_fract_data'] = [(round(p,2),[round(x,2) for x in ci]) for p, ci in zip(self.popgen_stats['p_poly_fract'], self.popgen_stats['poly_wilson'])]
                data_report['p_mono_unique'] = [(round(p,2),[round(x,2) for x in ci]) for p, ci in zip(self.popgen_stats['p_mono_unique'], self.popgen_stats['wilson_mono_unique'])]

                data_report['realmccoilcoi'] = [(round(p,2), [max(round(p - 1.96*std,2),0), round(p + 1.96*std,2)]) for p, std in zip(self.popgen_stats['mccoil_coi'], self.popgen_stats['mccoil_coi_std'])]
                data_report['realmccoilcoi_poly'] = [(round(p,2),[max(round(p - 1.96*std,2),0), round(p + 1.96*std,2)]) for p, std in zip(self.popgen_stats['mccoil_coi_poly'], self.popgen_stats['mccoil_coi_poly_std'])]

                data_report['RH_array'] = [(round(p,2),[round(p - 1.96*std,2), round(p + 1.96*std,2)]) for p, std in zip(self.RH_yearly_averages, np.sqrt(self.RH_yearly_variances))]
                data_report_df = pd.DataFrame.from_dict(data_report).T
                data_report_df.columns = self.chrono_years
                if output_file:
                    data_report_df.to_csv(output_file)
                return data_report_df

        file = "~{barcodes_file}"
        sheet_names = pd.ExcelFile(file).sheet_names
        regions = []
        for sheet in sheet_names:
            if sheet != "Control":  # Exclude the "Control" sheet
                df = pd.read_excel(file, sheet_name=sheet, nrows=1)  # Read only the first row for efficiency
                if 'mccoil_median' in df.columns:
                    regions.append(sheet)
                else:
                    continue
        print(regions)


        #-------------------------------------------------------------
        #Read in the data 
        BS = {}
        frames = []
        for ISO3 in regions:
            BS[ISO3.split(':')[0]] = Barcode_Stats(file, ISO3, 
                                    sheet_name= ISO3.replace(':','_'), adjustedN=False)#adjustedN=True)
            scramble_idxes = create_jackknife_lists(BS[ISO3.split(':')[0]].master_df)
            for k in range(1,10):
                print('prepping {ISO3} k_idx = {k}'.format(ISO3= ISO3, k = k))
                BS[ISO3.split(':')[0] + '.' + str(k)] = Barcode_Stats(file, ISO3, 
                                    sheet_name= ISO3.replace(':','_'), adjustedN=False, k_idx = k, scramble_idxes= scramble_idxes)#adjustedN=True)
            #frames.append(BS[ISO3.split(':')[0]].mono_barcode_df)
        #master_mono_df = pd.concat(frames)

        sample_distribution = []

        for ISO3 in BS.keys():
            sample_distribution += list(BS[ISO3].n_total)
            print(ISO3, BS[ISO3].n_total)
            

        #-------------------------------------------------------------
        #create covariate table

        def setup_covariate_table(k_idx = None):
            data_df = defaultdict(list)
            for ISO3 in [x.split(':')[0] for x in regions]:
                if k_idx:
                    ISO3 = ISO3 + '.' + str(k_idx)
                data_df['Region'] += list([ISO3 for _ in BS[ISO3].chrono_years])
                data_df['Year'] += list([year for year in BS[ISO3].chrono_years])
            # data_df['Incidence'] += list([incidences_dict[ISO3][year] for year in BS[ISO3].chrono_years])
                data_df['RH'] += list(BS[ISO3].RH_yearly_averages)
                data_df['RH_var'] += list(BS[ISO3].RH_yearly_variances)

                data_df['cotx'] += list(BS[ISO3].popgen_stats['cotx' + '_average'])
                data_df['cotx_var'] += list(BS[ISO3].popgen_stats['cotx_var'])

                data_df['monoclonality'] += list(BS[ISO3].popgen_stats['p_mono_clonal'])
                data_df['monoclonality_var'] += list(BS[ISO3].popgen_stats['var_mono_clonal'])
                data_df['monoclonality_ci'] += list(BS[ISO3].popgen_stats['wilson_mono_clonal'])

                data_df['poly_fract'] += list(BS[ISO3].popgen_stats['p_poly_fract'])
                data_df['poly_fract_var'] += list(BS[ISO3].popgen_stats['var_poly_fract'])
                data_df['poly_fract_ci'] += list(BS[ISO3].popgen_stats['poly_wilson'])

                data_df['mccoil_coi'] += list(BS[ISO3].popgen_stats['mccoil_coi'])
                data_df['mccoil_coi_std'] += list(BS[ISO3].popgen_stats['mccoil_coi_std'])
                data_df['mccoil_coi_var'] += list(np.asarray(BS[ISO3].popgen_stats['mccoil_coi_std'])**2)

                data_df['mccoil_coi_poly'] += list(BS[ISO3].popgen_stats['mccoil_coi_poly'])
                data_df['n'] += list(BS[ISO3].n_total)
                data_df['is_cachement'] += list([False for _ in BS[ISO3].chrono_years])

            # data_df['monoclonality_downsampled'] += list([x[0] for x in prop_clonal_arrays[ISO3]])
                #data_df['monoclonality_downsampled_ci'] += list([x[1] for x in prop_clonal_arrays[ISO3]])
            return DataFrame(data_df)

        data_dfs = {}
        data_dfs['all'] = setup_covariate_table()
        for k in range(0,10):
            data_dfs[k] = setup_covariate_table(k)
            
        #create_jackknife_estimates
        covariates = ["poly_fract", "monoclonality", "RH", "cotx", "mccoil_coi"]
        jackknife_results = np.stack([data_dfs[k][['Region', 'Year'] + covariates].to_numpy() for k in range(1,10)], axis = 2)
        jacknife_estimates = []
        for site_year_results in jackknife_results:
            site = site_year_results[0][0].split('.')[0]
            jackknife_estimate = [site] + list(np.mean(site_year_results[1:], axis = 1))
            jacknife_estimates.append(jackknife_estimate)
        jacknife_estimate_df = DataFrame(jacknife_estimates, columns = ['Region', 'Year'] + covariates)
        jacknife_estimate_df['is_cachement'] = False
        #-------------------------------------------------------------
        #run and combine predictions on all k-fold subsets for all leave-one-out population models
        with open("/usr/local/bin/regression_models.pkl", "rb") as file:  # "rb" stands for read-binary
            model = dill.load(file)

        #with open('retrain_model_coi_poly_mono_only.pkl', "rb") as file:  # "rb" stands for read-binary
        #    model = dill.load(file)


        def setup_model_predictions(subset_df):
            x_of_interest = ["Year", "poly_fract", "monoclonality", "RH", "cotx", "mccoil_coi", "is_cachement"]
            regions = list(subset_df['Region'])
            subset_df = subset_df[x_of_interest].copy()
            subset_df = pd.DataFrame(subset_df)
            subset_df['Incidence'] = 0
            subset_df['Region'] = regions
            subset_df.insert(0, 'Region', subset_df.pop('Region'))  # Ensure 'region' is the first column

            # Mapping model index to label
            model_mapping = {-10: "low", 10: "high"}
            data_M = {}
            # Iterate through both models (-10 and 10)
            for model_idx, label in model_mapping.items():
                response, predictor = dmatrices(model[model_idx].formula, subset_df, return_type="dataframe")
                X = np.asarray(predictor)

                # Get predictions from the main GLM model
                subset_df[f'Incidence_{label}'] = model[model_idx].glm_model.predict(X)

                # Store LOO model predictions
                model_data = {x: model[model_idx].loo_models[x].predict(X) for x in model[model_idx].loo_models.keys()}

                # Compute mean and confidence intervals
                mean_values, lower_ci, upper_ci = [], [], []
                model_predictions = []

                data_M[label] = np.asarray(list(model_data.values())).T
            return data_M

        #------------------------------------------------------------- 
        #collect results
        high_combined_pred = setup_model_predictions(jacknife_estimate_df)['high']
        low_combined_pred = setup_model_predictions(jacknife_estimate_df)['low']

        report_df = jacknife_estimate_df[["Region", "Year", "poly_fract", "monoclonality", "RH", "cotx", "mccoil_coi", "is_cachement"]]
        report_df['Incidence_high'] = np.mean(high_combined_pred, axis = 1)
        high_ci_upper = np.percentile(high_combined_pred, 97.5, axis = 1)
        high_ci_lower =  np.percentile(high_combined_pred, 2.5, axis = 1)
        report_df['Incidence_high_ci'] = [(round(x,4), round(y,4)) for x,y in zip(high_ci_lower, high_ci_upper)]

        report_df['Incidence_low'] = np.mean(low_combined_pred, axis = 1)
        low_ci_upper = np.percentile(low_combined_pred, 97.5, axis = 1)
        low_ci_lower =  np.percentile(low_combined_pred, 2.5, axis = 1)
        report_df['Incidence_low_ci'] = [(round(x,4), round(y,4)) for x,y in zip(low_ci_lower, low_ci_upper)]

        report_df = report_df.rename(columns={'monoclonality': 'monogenomic_clone_fraction'})

        model_data = {}
        for pop, sim_high, sim_low in zip(report_df['Region'], high_combined_pred, low_combined_pred):
            model_data[pop] = (list(sim_high), list(sim_low))
        # Save results to CSV
        report_df.to_csv("regression_output.csv", index=False)

        with open('raw_model_results.json','w') as fp:
            json.dump(model_data, fp)

        with open('data_dfs.pkl', 'wb') as fp:
            dill.dump(data_dfs,fp)

        EOF
    >>>

    output {
        File regression_output = "regression_output.csv"
        File covariate_tables = "data_dfs.pkl"
        File raw_model_results = "raw_model_results.json"
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/interpret-barcodes-python:v5"  
    }
}

workflow RunRegression {
    input {
        File barcodes_file
    }

    call RegressionOutput {
        input:
            barcodes_file = barcodes_file
    }

    output {
        File regression_output = RegressionOutput.regression_output
        File covariate_tables = RegressionOutput.covariate_tables
        File raw_model_results = RegressionOutput.raw_model_results
    }
}