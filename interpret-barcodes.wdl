version 1.0

task CompileBarcodeFiles {
    input {
        Array[File] input_directory
        Int threshold
    }

    Int disk_size = 1 + 5 * ceil(
    size(input_directory, "GB")
    )

    command <<<

        python3 <<EOF

import pandas as pd
from pandas import DataFrame, read_csv, read_excel
from collections import defaultdict, Counter
import re
import numpy as np
import glob
import copy
import sys
import os
import subprocess

assay_names = ['A1', 'B1', 'A2', 'B2', 'A3', 'B3', 'A4',
               'B4', 'A5', 'B5', 'A6', 'B6', 'A7', 'B7',
               'A8', 'B8', 'A9', 'B9', 'A10', 'B10', 'A11',
                'B11', 'A12', 'B12']
    

def bootstrap(array, iterations):
    resample = np.random.choice(array, (iterations, len(array)), replace=True)
    return np.mean(resample, axis = 1)

class Barcode_Sample:
    assay_names = ['A1', 'B1', 'A2', 'B2', 'A3', 'B3', 'A4',
               'B4', 'A5', 'B5', 'A6', 'B6', 'A7', 'B7',
               'A8', 'B8', 'A9', 'B9', 'A10', 'B10', 'A11',
                'B11', 'A12', 'B12']
    
    def __init__(self, raw_name, crt_df, setup_df, duplicate = False, 
                 a1_crt_field = 'Allele1 Crt',a2_crt_field = 'Allele2 Crt', full_sweep = False):
        if duplicate:
            raw_name, iteration = raw_name.split('_v')
        else:
            iteration = 1
        self.raw_name = raw_name
        self.cleaned_name = None
        self.wgs_name = None
            
    
        self.crt_df = crt_df
        self.assay_df = setup_df
        self.iteration = iteration
        self.reconstructed_name = None
        self.a1_crt_field = a1_crt_field
        self.a2_crt_field = a2_crt_field
        
        
        
        self.process_crt_df()
        
        self.barcode, self.delta_ct = {},[]
        
        for threshold in range(0,20):
            self.call_barcode(threshold_crt = threshold)
        if full_sweep:
            for threshold in np.linspace(0,20,100):
                self.call_barcode(threshold_crt = threshold)        
        self.return_n_poly()
        self.return_n_missing()
        self.call_mono_vs_poly()
    
    #def __repr__(self):
    #    return self.__dict__.keys()
    
    
    def process_crt_df(self):
        self.crt_df[self.a1_crt_field] = [float(x) if x != 'Undetermined' else 0 for x in self.crt_df[self.a1_crt_field].to_numpy()]
        self.crt_df[self.a2_crt_field] = [float(x) if x != 'Undetermined' else 0 for x in self.crt_df[self.a2_crt_field].to_numpy()]

        self.crt_df['Allele_Crt_Ratio'] = self.crt_df[self.a1_crt_field].to_numpy() / self.crt_df[self.a2_crt_field].to_numpy()
        self.crt_df['Allele_Crt_diff'] = np.absolute(self.crt_df[self.a1_crt_field].to_numpy() - self.crt_df[self.a2_crt_field].to_numpy())
    
    def call_barcode(self, min_crt = 0, threshold_crt = 5, maximum_crt = 38):
        barcode = {}
        for row in self.crt_df[['SNP Assay Name', self.a1_crt_field, self.a2_crt_field, 
                                   'Allele_Crt_Ratio']].to_numpy():
            assay_name, a1_crt, a2_crt, allele_ratio = row
            
            if (a1_crt == min_crt) & (a2_crt == min_crt):
                    call = 'X'
            elif allele_ratio == np.inf:
                #a1 call
                call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele1 Name'].to_list()[0]
            elif allele_ratio == 0:
                call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele2 Name'].to_list()[0]
            elif (a1_crt < maximum_crt) and (a2_crt < maximum_crt):
                if np.absolute(a1_crt - a2_crt) <= threshold_crt: 
                    call = 'N'

                elif a1_crt < a2_crt:
                    call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele1 Name'].to_list()[0]
                elif a1_crt > a2_crt:
                    call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele2 Name'].to_list()[0]
                
                self.delta_ct.append(np.absolute(a1_crt - a2_crt))

            elif a1_crt > maximum_crt:
                call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele2 Name'].to_list()[0]
                self.delta_ct.append(np.absolute(a1_crt - a2_crt))
            elif a2_crt > maximum_crt:
                call = self.assay_df[self.assay_df['SNP Assay Name'] == assay_name]['Allele1 Name'].to_list()[0]
                self.delta_ct.append(np.absolute(a1_crt - a2_crt))
            else:
                call = '?'
                

            barcode[assay_name] = call
                        
        self.barcode[threshold_crt] = np.asarray([barcode[k] if k in barcode else '?' for k in self.assay_names])
    
    def return_n_poly(self):
        self.poly_count = {}
        for threshold in self.barcode:
            counts = Counter(self.barcode[threshold])
            self.poly_count[threshold] = counts.pop('N', 0)
            
    def return_n_missing(self):
        self.n_missing = {}
        for threshold in self.barcode:
            counts = Counter(self.barcode[threshold])
            self.n_missing[threshold]=  counts.pop('X', 0)
    def call_mono_vs_poly(self, poly_threshold = 2):
        self.status = {}
        for threshold in self.poly_count:
            if self.poly_count[threshold] >= poly_threshold:
                self.status[threshold] = 'P'
            else:
                self.status[threshold] = 'M'
    
    
def extract_barcodes(file, barcode_dict=None, duplicates = defaultdict(lambda:1), output_control = False,
                    a1_crt_field = 'Allele1 Crt', a2_crt_field = 'Allele2 Crt', full_sweep = False):
    fin = read_excel(file, sheet_name = ['Results', 'Sample Setup'])
    header_idx = fin['Results'][fin['Results']['Block Type'] == 'Well'].index[0]
    results_header = fin['Results'].iloc[header_idx]
    setup_header = fin['Sample Setup'].iloc[header_idx]

    for sheet_name in ['Results', 'Sample Setup']:
        header = fin[sheet_name].iloc[header_idx]
        fin[sheet_name] = fin[sheet_name].iloc[header_idx + 1:]
        fin[sheet_name].columns = header
        fin[sheet_name].reset_index(inplace=True, drop = True)
    
    sample_setup_groups = {k:v for k,v in fin['Sample Setup'].groupby('Sample Name')}
    results_groups = {k:v for k,v in fin['Results'].groupby('Sample Name')}
    
    if not barcode_dict:
        barcode_dict = {}
        
    if not output_control:
        for raw_name in results_groups:
            if ('Row' not in raw_name):
                if raw_name in barcode_dict:
                    duplicates[raw_name] += 1
                    dup_name = raw_name + '_v{i}'.format(i= duplicates[raw_name])
                    barcode_dict[dup_name] = Barcode_Sample(dup_name, 
                                                                 results_groups[raw_name], 
                                                                 sample_setup_groups[raw_name], duplicate=True,
                                                           a1_crt_field= a1_crt_field, a2_crt_field = a2_crt_field, full_sweep = full_sweep)
                else:
                    barcode_dict[raw_name] = Barcode_Sample(raw_name, 
                                                                 results_groups[raw_name], 
                                                                 sample_setup_groups[raw_name],
                                                           a1_crt_field= a1_crt_field, a2_crt_field = a2_crt_field, full_sweep = full_sweep)
    
    else:
        for raw_name in results_groups:
            if raw_name in []:
                if raw_name in barcode_dict:
                    duplicates[raw_name] += 1
                    dup_name = raw_name + '_v{i}'.format(i= duplicates[raw_name])
                    barcode_dict[dup_name] = Barcode_Sample(dup_name, 
                                                                 results_groups[raw_name], 
                                                                 sample_setup_groups[raw_name], duplicate=True, full_sweep = full_sweep)
                else:
                    barcode_dict[raw_name] = Barcode_Sample(raw_name, 
                                                                 results_groups[raw_name], 
                                                                 sample_setup_groups[raw_name], full_sweep = full_sweep)
    return barcode_dict, duplicates



def read_barcode_files(directory=None, exclude_files= None):
    
    if not exclude_files:
        exclude_files = []
    tmp_dict = {}
    duplicates = defaultdict(lambda:1)
    if directory:
        pathname = directory + '/*xls*'
    else:
        pathname = '*xls*'
    for file in glob.glob(pathname):
        if file != 'barcodes.xlsx':
            print(file)
            if file not in exclude_files:
                try:
                    tmp_dict, duplicates = extract_barcodes(file, barcode_dict = tmp_dict, duplicates = duplicates)
                except:
                    tmp_dict, duplicates = extract_barcodes(file, barcode_dict = tmp_dict, duplicates = duplicates,
                                                                a1_crt_field='Allele1 Ct', a2_crt_field = 'Allele2 Ct')
    return tmp_dict
        
def prepare_output(sample_barcode_dict, remove_sample = [], threshold = 8):
    output_data = [['cc', 'ISO3', 'Year', 'Number_Text', 'Sample_Name', 'Raw_Name',
             'Barcode_String'] + assay_names + ['X', 'N', 'M/P', 'Delta_CT_Threshold']]
    
    current_path = os.getcwd()
    parent_path = os.path.basename(current_path)

    for sample in sample_barcode_dict.values():
        if sample not in remove_sample:
            if sample.raw_name == 'DD2':
                sample_name = sample.raw_name
                raw_name = sample.raw_name
                cc = ''
                year = parent_path
                iso3 = ''
                number_text = ''
            else:
                # Match patterns like 'KDG13' or 'SLP_20_7'
                match = re.match(r'([A-Z]+)(\d+)$|([A-Z]+)_(?:\d+_)?(\d+)', sample.raw_name)
                if match:
                    if match.group(1):  
                        iso3 = match.group(1)
                        number = match.group(2)
                    elif match.group(3): 
                        iso3 = match.group(3)
                        number = match.group(4)

                    number_text = str(number).zfill(3)
                    year = parent_path
                    sample_name = f"SEN_{year}_{iso3}_{number_text}"
                    raw_name = sample.raw_name
                    cc = 'SEN'
                else:
                    sample_name = sample.raw_name
                    raw_name = sample.raw_name
                    cc = ''
                    year = parent_path
                    iso3 = ''
                    number_text = ''

            line = [cc, iso3, year, number_text, sample_name, raw_name, ''.join(sample.barcode[threshold])] + list(sample.barcode[threshold]) + [sample.n_missing[threshold], sample.poly_count[threshold], sample.status[threshold], threshold]
            output_data.append(line)

    duplicate_counter = Counter([x[4] for x in output_data])
    duplicates = [c for c in duplicate_counter if duplicate_counter[c] != 1]
    final_df =DataFrame(output_data[1:], columns = output_data[0])  
    drop_indexes=[]
    for duplicate in duplicates:
        if duplicate and ('3D7' not in duplicate) and ('Dd2' not in duplicate):
            tmp_df = final_df[final_df['Sample_Name'] == duplicate]
            keep_index = tmp_df.index[np.argmin(tmp_df['X'])]
            drop_indexes +=[index for index in tmp_df.index if index != keep_index]
    final_df.drop(drop_indexes, inplace=True)
    return final_df

if __name__ == '__main__':
    barcode_data = read_barcode_files(directory="~{input_directory}")
    combined_df = prepare_output(barcode_data, threshold=~{threshold})
    combined_df.to_excel("compiled_output.xlsx", index=False)
EOF
    >>>
    
    output {
        File output_file = "compiled_output.xlsx"
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "python:3.8-slim"  
    }
}

workflow compile_datasets {
    input {
        Array[File] input_directory   
        Int threshold = 8       # Threshold value (default: 8)
    }

    call CompileBarcodeFiles {
        input:
            input_directory = input_directory,
            threshold = threshold
    }

    output {
        File combined_file = CompileBarcodeFiles.output_file
    }
}