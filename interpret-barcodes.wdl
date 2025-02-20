version 1.0

task CompileBarcodeFiles {
    input {
        Array[File] input_directory
        Int year
        Int threshold
    }

    Int disk_size = 1 + 5 * ceil(
    size(input_directory, "GB")
    )

    command <<<

        python3 <<EOF

        import pandas as pd
        from pandas import DataFrame, read_csv, read_excel
        import collections
        from collections import defaultdict, Counter
        import re
        import xlsxwriter
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



        def read_barcode_files(file_list=None, exclude_files=None):
            if not exclude_files:
                exclude_files = []
            if not file_list:
                raise ValueError("No input files provided. Please supply a list of Excel files.")

            tmp_dict = {}
            duplicates = defaultdict(lambda: 1)

            for file in file_list:
                if file != 'barcodes.xlsx' and file not in exclude_files:
                    print(f"Processing file: {file}")
                    try:
                        tmp_dict, duplicates = extract_barcodes(
                            file, barcode_dict=tmp_dict, duplicates=duplicates
                        )
                    except Exception as e:
                        print(f"Error processing {file}: {e}. Attempting with default fields.")
                        tmp_dict, duplicates = extract_barcodes(
                            file,
                            barcode_dict=tmp_dict,
                            duplicates=duplicates,
                            a1_crt_field='Allele1 Ct',
                            a2_crt_field='Allele2 Ct'
                        )

            return tmp_dict

        def prepare_output(sample_barcode_dict, remove_sample = [], threshold = 8):
            output_data = [['cc', 'ISO3', 'Year', 'Number_Text', 'Sample_Name', 'Raw_Name',
                    'Barcode_String'] + assay_names + ['X', 'N', 'M/P', 'Delta_CT_Threshold']]
            

            for sample in sample_barcode_dict.values():
                if sample not in remove_sample:
                    raw_name = sample.raw_name

                    if raw_name == 'DD2':
                        # Handle Dd2
                        sample_name = raw_name
                        raw_name = raw_name
                        cc = ''
                        year = ~{year}
                        iso3 = ''
                        number_text = ''

                    else:
                        # Match 'SEN_2023_KDG.P_018' format
                        if match := re.match(r'^(SEN)_(\d{4})_([A-Z]+)(\.[A-Z]+)?_(\d+)$', raw_name):
                            cc, year, iso3, suffix, number_text = match.groups()
                            number_text = number_text.zfill(3)
                            #suffix = suffix or ''   Ensure suffix is not None
                            sample_name = raw_name

                        # Match 'KDG018', 'BAN113', or 'BAN113.P'
                        elif match := re.match(r'([A-Z]+)(\d+)(\.[A-Z]+)?$', raw_name):
                            iso3 = match.group(1)  # Extract 'KDG' or 'BAN'
                            number = match.group(2)  # Extract '018' or '113'
                            suffix = match.group(3) or ''  # Extract '.P' or empty if not present
                            number_text = str(number).zfill(3)
                            year = ~{year}
                            sample_name = f"SEN_{year}_{iso3}_{number_text}"  # Keep original formatting
                            cc = 'SEN'

                        # Match 'SLP_20_7' format
                        elif match := re.match(r'([A-Z]+)_(?:\d+_)?(\d+)(\.[A-Z]+)?$', raw_name):
                            iso3 = match.group(1)  # Extract 'SLP'
                            number = match.group(2)  # Extract '7'
                            suffix = match.group(3) or ''  # Extract '.P' or empty if not present
                            number_text = str(number).zfill(3)
                            year = ~{year}
                            sample_name = f"SEN_{year}_{iso3}_{number_text}"
                            cc = 'SEN'

                        else:
                            # Fallback logic for unknown formats
                            sample_name = raw_name
                            raw_name = raw_name
                            cc = ''
                            year = ~{year}
                            iso3 = ''
                            number_text = ''

                    # Build the line and append to output_data
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
            barcode_data = read_barcode_files(file_list = ["~{sep='","' input_directory}"])
            out = prepare_output(barcode_data, threshold=~{threshold})
            with pd.ExcelWriter('barcodes.xlsx', engine='xlsxwriter') as writer:
                for iso3, group in out.groupby('ISO3'):
                    if pd.isna(iso3) or iso3 == '':
                        sheet_name = 'Control'
                    else:
                        sheet_name = str(iso3)
                    
                    group.to_excel(writer, sheet_name=sheet_name, index=False)
    
        EOF
    >>>
    
    output {
        File output_file1 = 'barcodes.xlsx'
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/interpret-barcodes-python:v3"  
    }
}

task CreateMcCoilExcel {
    input {
        File output_file1
    }

    Int disk_size = 1 + 5 * ceil(
    size(output_file1, "GB")
    )

    command <<<

        python3 <<EOF

        import pandas as pd
        from pandas import DataFrame, read_csv, read_excel
        import collections
        from collections import defaultdict, Counter
        import re
        import xlsxwriter
        import numpy as np
        import copy
        import sys
        import os
        import subprocess

        recoded_barcodes_dict = {}
        usable_sites = []
        with pd.ExcelFile("~{output_file1}") as xcel_file:
            for sheetname in xcel_file.sheet_names:
                barcodes_df = xcel_file.parse(sheetname)
                pop_id = sheetname
                if len(barcodes_df) > 1: # Should be greater than 30, argument as script, will make this a parameter in Terra!
                    #print(pop_id, len(barcodes_df))
                    usable_sites.append(pop_id)
                    barcodes_df = barcodes_df[(barcodes_df['ISO3'].isin([pop_id]))]

                    start_idx = 7
                    end_idx = 31

                    # Process the loci columns
                    loci_position_names = barcodes_df.columns[start_idx:end_idx]
                    for position in loci_position_names:
                        barcodes_df[position] = [x.strip().upper() for x in barcodes_df[position]]  # Fix tab formatting
                    barcodes_df['M/P'] = [x.strip() for x in barcodes_df['M/P']]

                    # Filter mono barcode dataframe
                    mono_barcode_df = barcodes_df[(barcodes_df['X'] <= 2) & (barcodes_df['M/P'] == 'M')]
                    chrono_years = sorted(Counter(barcodes_df['Year']).keys())
                    loci_allele_dict = {}

                    # Get major and minor alleles
                    for column in mono_barcode_df.columns[start_idx:end_idx]:
                        counts = Counter(mono_barcode_df[column].to_numpy())
                        counts.pop('X', 0)
                        counts.pop('N', 0)
                        if counts:
                            major_allele = max(counts, key=counts.get)
                            counts.pop(major_allele, 0)
                            try:
                                minor_allele = max(counts, key=counts.get)
                            except:
                                minor_allele = None
                            loci_allele_dict[column] = major_allele
                        else:
                            # Handle empty counts case (e.g., set major_allele to None or skip this column)
                            loci_allele_dict[column] = None

                    # Recode the barcodes
                    loci_columns = list(barcodes_df.columns[start_idx:end_idx])
                    all_barcodes = barcodes_df[['Sample_Name'] + loci_columns].to_numpy()

                    recoded_barcodes = [['ind_name'] + list(loci_allele_dict.keys())]
                    for row in all_barcodes:
                        samplename = row[0]
                        barcode = row[1:]
                        recoded_barcode = [samplename]
                        for position, assay in zip(loci_allele_dict.keys(), barcode):
                            major = loci_allele_dict[position]
                            if major == assay:
                                recode = 1
                            elif assay == 'N':
                                recode = 0.5
                            elif assay == 'X':
                                recode = -1
                            else:
                                recode = 0
                            recoded_barcode.append(recode)
                        recoded_barcodes.append(recoded_barcode)
                    
                    recoded_barcodes_dict[pop_id] = pd.DataFrame(np.asarray(recoded_barcodes)[1:],
                                                                columns=np.asarray(recoded_barcodes)[0])

        with pd.ExcelWriter('barcodes_mccoil.xlsx', engine='xlsxwriter') as writer:
            for pop_id in recoded_barcodes_dict:
                recoded_barcodes_dict[pop_id].to_excel(writer, 
                                                    sheet_name=str(pop_id).replace(':', '_'), 
                                                    index=False, header=True)
        EOF
    >>>
    
    output {
        File output_file2 = 'barcodes_mccoil.xlsx'
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/interpret-barcodes-python:v3"  
    }
}

task RunMcCoil {
    input {
        File output_file2
    }

    Int disk_size = 1 + 5 * ceil(
    size(output_file2, "GB")
    )

    command <<<
        Rscript -e '
        library(readxl);
        execution_dir <- getwd();
        setwd("/THEREALMcCOIL/categorical_method");
        source("McCOIL_categorical.R");
        print("done");
        if (file.exists("McCOIL_categorical_code.so")) {
            dyn.load("McCOIL_categorical_code.so");
        } else {
            stop("Shared object file not found.");
        }
        excel_file = "~{output_file2}";
        sheet_names <- excel_sheets(excel_file);
        sheets_to_process <- setdiff(sheet_names, "Control");
        output_files <- c();
        output_summary_files <- c();
        for (sheet in sheets_to_process) {
            barcodes_mccoil <- read_excel(excel_file, sheet = sheet);
            if (nrow(barcodes_mccoil) < 20) {
                cat(paste("Skipping sheet", sheet, "because it has fewer than 20 observations."));
                next;
            }
            data <- as.data.frame(barcodes_mccoil[,-1]);
            rownames(data) <- barcodes_mccoil$ind_name;
            name <- paste0("mccoil_", sheet);

            McCOIL_categorical(
                data = as.matrix(data),
                maxCOI = 25,
                threshold_ind = 20,
                threshold_site = 5, 
                # SHOULD BE 20 for threshold_site
                totalrun = 1000,
                burnin = 100,
                M0 = 15,
                e1 = 0.05,
                e2 = 0.05,
                err_method = 3,
                path = getwd(),
                output = name);
    
            output_file <- name;
            output_summary_file <- paste0(name, "_summary.txt");
            file.copy(file.path(getwd(), output_file), file.path(execution_dir, output_file), overwrite = TRUE);
            file.copy(file.path(getwd(), output_summary_file), file.path(execution_dir, output_summary_file), overwrite = TRUE)
            }
        ' 
    >>>
        
    output {
        Array[File] mccoil_output_files = glob("mccoil_*")
        Array[File] mccoil_summary_files = glob("mccoil_*_summary.txt")
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/mccoil-docker:v5"  
    }
}

task FinalDataset {
    input {
        File output_file1
        Array[File] mccoil_summary_files
    }

    Int disk_size = 1 + 5 * ceil(
    size(output_file1, "GB")
    )

    command <<<

        python3 <<EOF

        import pandas as pd
        from pandas import DataFrame, read_csv, read_excel
        import collections
        from collections import defaultdict, Counter
        import re
        import xlsxwriter
        import numpy as np
        import copy
        import sys
        import os
        import subprocess
        import pathlib as Path

        mccoil_summaries= ["~{sep='","' mccoil_summary_files}"]

        with pd.ExcelFile("~{output_file1}") as xls:
            all_sheets = {sheet: xls.parse(sheet) for sheet in xls.sheet_names}

        good_sheets = []

        sheet_names = all_sheets.keys()

        for x in sheet_names:
            if x == 'Control':
                continue
            if len(all_sheets[x]) >= 20:
                good_sheets.append(x)

        for sheet_name in good_sheets:

            txt_file = next((f for f in mccoil_summaries if f"mccoil_{sheet_name}_summary.txt" in f), None)

            txt_data = pd.read_csv(txt_file, sep="\t")
            
            # Rename 'name' in txt_data to match 'Sample_Name' in Excel
            if 'name' in txt_data.columns:
                txt_data.rename(columns={'name': 'Sample_Name'}, inplace=True)

            # Merge on Sample_Name
            merged_data = all_sheets[sheet_name].merge(txt_data[['Sample_Name', 'median']], 
                                                    on='Sample_Name', how='left')

            # Rename the median column to 'mccoil_median'
            merged_data.rename(columns={'median': 'mccoil_median'}, inplace=True)

            # Update the dictionary with merged data
            all_sheets[sheet_name] = merged_data

        with pd.ExcelWriter('barcode_final.xlsx', engine='xlsxwriter') as writer:
            for sheet_name, sheet_data in all_sheets.items():
                sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
        EOF
    >>>

    output {
        File final_output = 'barcode_final.xlsx'
    }

    runtime {
        cpu: 2
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " SSD"
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 0
        docker: "karinab2000/interpret-barcodes-python:v3"  
    }
}

workflow compile_datasets {
    input {
        Array[File] input_directory   
        Int year 
        Int threshold = 8     
    }

    call CompileBarcodeFiles {
        input:
            input_directory = input_directory,
            year = year,
            threshold = threshold
    }

    call CreateMcCoilExcel {
        input:
            output_file1 = CompileBarcodeFiles.output_file1
    }

    call RunMcCoil {
        input:
            output_file2 = CreateMcCoilExcel.output_file2
    }

    call FinalDataset {
        input:
            output_file1 = CompileBarcodeFiles.output_file1,
            mccoil_summary_files = RunMcCoil.mccoil_summary_files
    }

    output {
        File final_output = FinalDataset.final_output
    }
}