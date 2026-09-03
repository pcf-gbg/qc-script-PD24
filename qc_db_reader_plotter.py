import datetime
import integrationutils as ius
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import pandas as pd
import seaborn as sns
sns.set()
import sqlite3
import sys
import warnings


''' Direct questions and concerns regarding this script to Egor Vorontsov
    egor.vorontsov@gu.se
'''

def bin_number(inlist,divider):
    span = max(inlist) - min(inlist)
##    print('The span of the mass scale is', span)
    return int(span/divider)

def create_reporting_fnames(qc_db_path, num_res, inst_name, raw_fname):
    ''' Creates the filenames for the graphical report and Excel output
        Based on the relative path which depends on the QC database path qc_db_path
        num_res: the number of the latest results fetched from the database
        inst_name: instrument name
        raw_fname: raw file name
        returns the dictionary:
        {'excel_out':full path, 'graph_out':full_path}
    '''
    qc_db_pobj = pathlib.Path(qc_db_path)
    excel_name = inst_name + '_latest_' + str(num_res) + '_results.xlsx'
    excel_path = qc_db_pobj.parents[1].joinpath(excel_name)
    png_name = 'qc_plot_' + raw_fname[:-4] + '.png'
    png_path = qc_db_pobj.parents[1].joinpath(png_name)

    return {'excel_out':excel_path, 'graph_out':png_path}

def find_hela_peptides(in_dfs):
    ''' Finds the retention times of the 7 selected peptides from the Peptide Groups table.
        Returns 45.0 min if the peptide is not found
    '''
    selected_hela_peptides = {'STELLIR':'pept_416',
                              'HLQLAIR':'pept_425',
                              'AGFAGDDAPR':'pept_488',
                              'VFLENVIR':'pept_495',
                              'DAVTYTEHAK':'pept_567',
                              'VAPEEHPVLLTEAPLNPK':'pept_652',
                              'TVTAMDVVYALK':'pept_655'
                              }
    try:
        # Select the Peptide Groups df from the dictionary
        df = in_dfs['Peptide Groups']
        # Find the column that contains the retention time
        for c in df.columns:
            if 'RT [min] (by Search Engine)' in c:
                rt_col = c
        rt_dict = {}
        for pep_seq in selected_hela_peptides.keys():
            df_pept = df[(df['Sequence'] == pep_seq) & (df['Modifications'].isnull())]
            if len(df_pept.index) == 1:
                rt_dict[selected_hela_peptides[pep_seq]] = df_pept[rt_col].iloc[0]
            else:
                rt_dict[selected_hela_peptides[pep_seq]] = None
        
        return rt_dict
    except:
        warnings.warn('Could not process Peptide Groups table for the extraction of the selected peptides',
                      UserWarning)
        return None

def fname_to_instrument(infname):
    '''Decides on the instrument name based on the file name'''
    infname = infname.lower()
    instrument_name = None
    if ('lumos_' in infname) and ('faims' in infname):
        instrument_name = 'Lumos_FAIMS'
    elif ('lumos_' in infname) and ('faims' not in infname):
        instrument_name = 'Lumos'
    elif ('fusion_' in infname) and ('faims' in infname):
        instrument_name = 'Fusion_FAIMS'
    elif ('fusion_' in infname) and ('faims' not in infname):
        instrument_name = 'Fusion'
    elif ('qehf_' in infname):
        instrument_name = 'QEHF'
    elif ('qe_' in infname):
        instrument_name = 'QE'
    elif ('elite_' in infname):
        instrument_name = 'Elite'
    elif ('exp_' in infname) and ('faims' not in infname):
        instrument_name = 'Exp'
    elif ('exp_' in infname) and ('faims' in infname):
        instrument_name = 'Exp_FAIMS'
    elif ('eclipse_' in infname) and ('faims' in infname):
        instrument_name = 'Eclipse_FAIMS'
    elif ('eclipse_' in infname) and ('faims' not in infname):
        instrument_name = 'Eclipse'

    return instrument_name

def generate_input_dict(in_dfs):
    ''' Reads the dictionary of the experimental tables.
        Generates the tuple for writing into the SQLite DB,
        As well as the arrays that are used for plotting.
        The output is a dictionary:
        {'Tuple for database':tuple,
        'PSM inj time array':np.array,
        'MSMS inj time array':np.array,
        'Prec int array':np.array,
        'Peak width array':np.array,
        'Prec mz array':np.array,
        'ppm error array':np.array,
        'RT array':np.array,
        'Search engine score array':np.array
        }
    '''
    rawfnames, rawfdates, instr_from_fn, _ = parse_table_input_file(in_dfs)
    search_time = datetime.datetime.now().isoformat(timespec='seconds')
    
    msms_inj_times,prec_ints_msms = parse_qc_table_msms(in_dfs)
    msms_num = len(msms_inj_times)
    mean_it = round(np.nanmean(msms_inj_times),3)
    median_it = round(np.nanmedian(msms_inj_times),3)
    
    (psm_mzs, psm_mhplus, psm_ppm_errs, psm_rts,
     se_scores, psm_inj_times, sengine, _) = parse_qc_table_psm(in_dfs)
    psm_num = len(psm_mzs)
    mean_psm_it = round(np.nanmean(psm_inj_times),3)
    median_psm_it = round(np.nanmedian(psm_inj_times),3)
    mean_error = round(np.nanmean(psm_ppm_errs),3)
    median_error = round(np.nanmedian(psm_ppm_errs),3)
    err_sdev = round(np.nanstd(psm_ppm_errs),3)
    mean_se_score = round(np.nanmean(se_scores),2)

    hela_pept_rts = find_hela_peptides(in_dfs)

    (mean_p_width, p_width_stdev, p_width_array,
     sum_feature_abund, lfq_int_array) = parse_cons_features(in_dfs)
    if lfq_int_array is None:
        prec_ints = prec_ints_msms
    else:
        prec_ints = lfq_int_array
    mean_prec_int = int(np.nanmean(prec_ints))
    
    pept_num = return_peptide_number(in_dfs)
    prot_num = return_protein_number(in_dfs)

    print('Input contained',prot_num,'Master proteins')
    print('Input contained',pept_num,'peptides')
    print('Input contained',psm_num,'PSMs')
    print('Input contained',msms_num,'MSMS spectra')
    
    out_tuple = (None, rawfnames[0], rawfdates[0], search_time, instr_from_fn[0],
                 prot_num, pept_num, psm_num, msms_num,
                 round((psm_num/msms_num),3),
                 mean_psm_it, median_psm_it, mean_it, median_it,
                 mean_error, median_error, err_sdev,
                 sum_feature_abund, mean_prec_int, mean_se_score,
                 mean_p_width, p_width_stdev,
                 hela_pept_rts['pept_416'],hela_pept_rts['pept_425'],hela_pept_rts['pept_488'],
                 hela_pept_rts['pept_495'],hela_pept_rts['pept_567'],hela_pept_rts['pept_652'],
                 hela_pept_rts['pept_655'], ('Searched with '+str(sengine)))

##    colnames = ('search_id','raw_file','file_date','search_date','instrument',
##                'protein_number','peptide_number','psm_number',
##                'msms_number','id_rate',
##                'mean_psm_it_ms','median_psm_it_ms',
##                'mean_msms_it_ms','median_msms_it_ms',
##                'mean_mz_err_ppm','median_mz_err_ppm','mz_err_ppm_stdev',
##                'total_prec_intensity','mean_prec_intensity',
##                'mean_sengine_score',
##                'mean_peak_width','peak_width_stdev',
##                'pept_416','pept_425','pept_488','pept_495',
##                'pept_567','pept_652','pept_655',
##                'comment')

    out_dict = {}
    out_dict['Tuple for database'] = out_tuple
    out_dict['PSM inj time array'] = psm_inj_times
    out_dict['MSMS inj time array'] = msms_inj_times
    out_dict['Prec int array'] = prec_ints
    out_dict['Peak width array'] = p_width_array
    out_dict['Prec mz array'] = psm_mzs
    out_dict['ppm error array'] = psm_ppm_errs
    out_dict['RT array'] = psm_rts
    out_dict['Search engine score array'] = se_scores
    
    return out_dict

def generate_simulated_tuple():

    fname = 'Nothing_' + str(np.random.randint(10,200)) + '.raw'
    d = '2019/10/' + str(np.random.randint(10,31))
    s = '2019/10/24'
    i = 'Eclipse'
    psmn = np.random.randint(8000,14000)
    pn = psmn - np.random.randint(50,500)
    msmsn = np.random.randint(25000,33000)
    mean_psm_it = 10 + 15*np.random.random()
    median_psm_it = mean_psm_it - np.random.randint(1,5)
    mean_it = 15 + 15*np.random.random()
    median_it = mean_it - np.random.randint(1,5)
    mean_er = np.random.randint(-5,5)*np.random.random()
    median_er = mean_er - 0.1*mean_er*np.random.random()
    er_dev = abs(mean_er)*0.5*np.random.random()
    mean_int = np.random.randint(600000,2000000)
    mean_score = 15 + 15*np.random.random()
    
    out_tuple = (None,fname,d,s,i,1800,pn,psmn,msmsn,psmn/msmsn,
                 mean_psm_it,median_psm_it,mean_it,median_it,
                 mean_er,median_er,er_dev,
                 mean_int,mean_score,
                 27,25,19.1,35.5,16.5,32.1,40.6,
                 'Simulated')
    return out_tuple

def get_console_arg():
    try:
        assert (len(sys.argv) == 2), "Script requires only one parameter"
        
        print('Input file',sys.argv[1])
        injs = sys.argv[1]
        print('Current work dir', os.getcwd())
        return injs
    except:
        print('Could not pass json file to the script')
        return None

def parse_cons_features(in_dfs):
    ''' If the Consensus Features table is present, returns the tuple with peak properties
        (for features with charge 2 or higher):
        (mean peak width in s, st dev of peak width in s,
        numpy array of peak widths in s, numpy array of feature intensity)
        If the Consensus Features table is absent, returns
        (None, None, numpy array of zeroes, None)
    '''
    if 'Consensus Features' in in_dfs.keys():
        try:
            df = in_dfs['Consensus Features']
            df1 = df[df['Charge'] > 1]
            diff_array = np.subtract(np.array(df1['Right RT [min]']),np.array(df1['Left RT [min]']))
            diff_array = 60*diff_array
            try:
                for c in df1.columns:
                    if 'Abundances (per File)' in c:
                        abundcol = c
                int_array = np.array(df1[abundcol])
            except:
                warnings.warn('Could not find Consensus Feature intinsities', UserWarning)
                int_array = None

            if int_array is None:
                sum_abund = None
            else:
                sum_abund = np.sum(int_array)
                print('Sum of all intensities is',sum_abund)
            
            return (round(np.nanmean(diff_array),1),
                    round(np.nanstd(diff_array),1),
                    diff_array,
                    sum_abund, int_array)
        except:
            warnings.warn('Could not process Consensus Features table', UserWarning)
            return (None, None, np.zeros(100,dtype=int), None, None)
    else:
        print('Could not find the Consensus Features table')
        return (None, None, np.zeros(100,dtype=int), None, None)

def parse_qc_table_msms(in_dfs):
    ''' Returns a tuple of numpy arrays
        (injection times array, precursor intens array)
        '''
    try:
        # Select the MSMS df from the dictionary
        df = in_dfs['MS/MS Spectrum Info']
        inj_times = np.array(df['Ion Inject Time [ms]'])
        prec_ints = np.array(df['Precursor Intensity'])

        return (inj_times,prec_ints)
    except:
        warnings.warn('Could not process MS/MS Spectrum Info table', UserWarning)
        return None

def parse_qc_table_psm(in_dfs):
    ''' Returns a tuple of numpy arrays
        (m/z values, MH+ mass values, ppm mass errors, ret times,
        search engine scores, search engine, ions matched per PSM)
        The search engine score to use dependes on what columns are present in the file.
        If neither Ions Score nor Xcorr column is found, returns an array made of zeros:
        sengine: 'Mascot' or 'Sequest', default is 'Mascot'
        '''

    try:
        # Select the score. Currently recognizes Mascot and Sequest
        df = in_dfs['PSMs']
        if 'Ions Score' in df.columns:
            seng_keyword = 'Ions Score'
            sengine = 'Mascot'
            print('Search engine identified as',sengine)
        elif 'XCorr' in df.columns:
            seng_keyword = 'XCorr'
            sengine = 'SequestHT'
            print('Search engine identified as',sengine)
        else:
            warnings.warn('Could not find the proper matching Score in the PSMs table', UserWarning)
            seng_keyword = None
            sengine = 'NA'
            print('Search engine NOT identified')
    except:
        warnings.warn('Could not find the proper matching Score in the PSMs table', UserWarning)
        seng_keyword = None
        sengine = 'NA'
        print('Search engine NOT identified')
    
    try:
        # Select the PSMs df from the dictionary
        df = in_dfs['PSMs']
        mzs = np.array(df['m/z [Da]'])
        mhplus = np.array(df['MH+ [Da]'])
        ppm_errors = np.array(df['DeltaM [ppm]'])
        rts = np.array(df['RT [min]'])
        inj_times = np.array(df['Ion Inject Time [ms]'])
        ions_matched = np.array(df['Matched Ions'])
        if seng_keyword is not None:
            se_scores = np.array(df[seng_keyword])
        elif seng_keyword is None:
            se_scores = np.zeros(len(mzs), dtype=int)

        return (mzs,mhplus,ppm_errors,rts,se_scores,
                inj_times,sengine,ions_matched)
    except:
        warnings.warn('Could not process PSMs table', UserWarning)
        return None

def parse_table_input_file(in_dfs):
    ''' Returns a tuple
        (instrument name list, creation date list, instr from name list, instr from metadata list)
        '''
    try:
        # Select the Input Files df from the dictionary
        df = in_dfs['Input Files']
        # Select the rows with .RAW files (df in principle can contain many .RAW files, .MSF etc)
        df1 = df[df['File Name'].str.contains('.raw')]
        # Create a list of base filenames for .RAW files
        shortnames = []
        for longname in df1['File Name']:  
            shortnames.append(pathlib.Path(longname).name)
        filedates = list(df1['Creation Date'])
        instr_from_metadata = list(df1['Instrument Name'])
        instr_from_fnames = []
        for n in shortnames:
            instr_from_fnames.append(fname_to_instrument(n))
        return (shortnames,filedates,instr_from_fnames,instr_from_metadata)
    except:
        warnings.warn('Could not process Input Files table', UserWarning)
        return None

def print_all_rows(conn, table):
    cur = conn.cursor()
    sql_command = "SELECT * FROM " + table
    cur.execute(sql_command)

    rows = cur.fetchall()
    for row in rows:
        print(row)
    cur.close()
    return True

def return_idrate_msms(in_df):
    '''Gets the whole QC dataframe
    Returns the long (melted) dataframe with the ID rate(%) and MS/MS in thousands'''
    new_df_dict = {}
    new_df_dict['search_id'] = np.array(in_df['search_id'])
    
    msms_array = np.array(in_df['msms_number'])
    msms_array = msms_array/1000
    new_df_dict['MSMS spectra/1000'] = msms_array

    id_rate_array = np.array(in_df['id_rate'])
    id_rate_array = id_rate_array*100
    new_df_dict['ID rate,%'] = id_rate_array
    
    out_df = pd.DataFrame(new_df_dict)
    
    out_df = pd.melt(out_df,id_vars=('search_id',),
                     value_vars=('MSMS spectra/1000','ID rate,%'),
                     var_name='Parameter',value_name='MSMS or IDrate')
    
    return out_df

def return_last_rows(conn, table, index_col, num, colnames):
    cur = conn.cursor()
    sql_command = ('''SELECT * FROM (
                        SELECT * FROM ''' + table + " ORDER BY " +
                   index_col + ''' DESC LIMIT ''' + str(num) + ''')
                    ORDER BY ''' + index_col + " ASC;")
    cur.execute(sql_command)

    rows = cur.fetchall()

    list_of_tuples = []
    for row in rows:
        list_of_tuples.append(row)

    cur.close()
    df = pd.DataFrame(list_of_tuples, columns=colnames)
    print('Fetched', num, 'latest results from the QC database')
    
    return df

def return_latest_psm_is(df, id_col, file_col, instr_col, psm_col):
    ''' Extracts info on PSM number, search ID and Instrument from the last row in DB
    '''
    last_row = df.iloc[-1]
    search_id = last_row[id_col]
    instr = last_row[instr_col]
    psm = last_row[psm_col]
    psm_string = str(psm) + ' PSMs in file ' + str(last_row[file_col])
    print('String to put on the graph', psm_string)
    return (search_id, instr, psm, psm_string)

def return_peptide_number(in_dfs):
    ''' Returns the number of peptides based on the Peptide table
        '''
    try:
        # Select the Peptide Groups df from the dictionary
        df = in_dfs['Peptide Groups']
        return (len(df.index))
    except:
        warnings.warn('Could not process Peptide Groups table', UserWarning)
        return None

def return_protein_number(in_dfs):
    ''' Returns the number of Master proteins based on the Proteins table
        '''
    try:
        # Select the Proteins df from the dictionary
        df = in_dfs['Proteins']
        # Select the Master proteins
        df1 = df[df['Master'] == 'IsMasterProtein']
        return (len(df1.index))
    except:
        warnings.warn('Could not process Proteins table', UserWarning)
        return None

def show_histo(inlist,labl,divider):
    binnum = bin_number(inlist,divider)
    ax = plt.hist(x=inlist, bins=binnum,color='blue',alpha=0.7,rwidth=0.85)
##    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(labl)
    plt.ylabel('Frequency')
##    plt.title('Wildcard Modification Masses from Byonic Search;'+
##              ' Searched with fixed TMT at N-termini and K, and fixed CAM on C')
##    plt.xticks(customticks)
##    plt.show()
    return True

def show_histo_fixedbinnum(inlist,labl,binnum):
    ax = plt.hist(x=inlist, bins=binnum,color='blue',alpha=0.7,rwidth=0.85)
    plt.xlabel(labl)
    plt.ylabel('Frequency')

    return True

def show_plots(df,fresh_data,plotsavingpath='D:\\plot.png',onlysave=False):
    ''' Shows the plots in the pre-defined layout.
        Requires the DataFrame df with thw latest N results from the database,
        And the fresh_data dictionary, that contains the necessary arrays from the new data.
        onlysave - bool, default is False
        If onlysave is False, the plot will be shown,
        If onlysave is True, the plot will be saved into the file plotsavingpath
    '''
    sid, latest_instr, psmnum, text_string = return_latest_psm_is(df, index_col, 'raw_file',
                                                                  'instrument', 'psm_number')
    f = plt.figure(figsize=(27,12),dpi=100)
    grid = plt.GridSpec(3, 7, figure=f, wspace=0.4, hspace=0.3)

    general_title = ('File ' + fresh_data['Tuple for database'][1] + 
                     ' from Instrument ' + latest_instr)
    plt.suptitle(general_title)

    # 1 Line plot of Peptide IDs and PSMs vs Search ID
    df_peptpsms = pd.melt(df,id_vars=('search_id',),
                              value_vars=('peptide_number','psm_number'),
                              var_name='Parameter',value_name='Peptides or PSMs')
    plt.subplot(grid[0, 0:2])
    ax1 = sns.lineplot(x='search_id', y='Peptides or PSMs', hue='Parameter',
                      data=df_peptpsms)
    ax1.legend(loc='lower left', fontsize='small', framealpha=0.5)
    
    plt.text(sid,psmnum,text_string,horizontalalignment='right',
             verticalalignment='bottom',
             fontdict={'color':'black','size':10})

    # 2 Line plot of the MSMS number and ID rate vs Search ID
    df_msmsid = return_idrate_msms(df)
    plt.subplot(grid[0, 2])
    ax2 = sns.lineplot(x='search_id', y='MSMS or IDrate', hue='Parameter',
                      data=df_msmsid)
    ax2.legend(loc='lower left', fontsize='small', framealpha=0.5)

    # 3 Line plot of the mean and median PSM injection time vs Search ID
    df_injtimes = pd.melt(df,id_vars=('search_id',),
                          value_vars=('mean_psm_it_ms', 'median_psm_it_ms'),
                          var_name='Parameter',value_name='Inj time, ms')
    plt.subplot(grid[0, 3])
    ax3 = sns.lineplot(x='search_id', y='Inj time, ms', hue='Parameter',
                      data=df_injtimes)
    ax3.legend(loc='lower left', fontsize='small', framealpha=0.5)

    # 4 Histogram of MSMS injection times for PSMs
    plt.subplot(grid[0, 4])
    # Use 5 to get bins of approximately 5 ms
    show_histo(fresh_data['PSM inj time array'],'PSM inj times, ms',5)
    
    # 5 Histogram of all the MSMS injection times
    plt.subplot(grid[0, 5])
    # Use 5 to get bins of approximately 5 ms
    show_histo(fresh_data['MSMS inj time array'],'All MSMS inj times, ms',5)

    # 6 Line plot of the Retention Times of the selected HeLa peptides vs Search ID
    df_hela_rts = pd.melt(df,id_vars=('search_id',),
                          value_vars=('pept_416','pept_425','pept_488',
                                      'pept_495','pept_567','pept_652','pept_655'),
                          var_name='Peptide',value_name='RT, min')
    plt.subplot(grid[0:2, 6])
    ax6 = sns.lineplot(x='search_id', y='RT, min', hue='Peptide',
                      data=df_hela_rts)
    ax6.legend(loc='lower left', fontsize='x-small', framealpha=0.5)
    
    # 7 Scatter plot of ppm mass error vs precursor m/z
    plt.subplot(grid[1:, 0])
    ax7 = sns.lineplot(fresh_data['Prec mz array'],fresh_data['ppm error array'],
                          alpha=0.4,s=10,c=('#1EA922',),marker='X')
    ax7.set_xlabel('Precursor m/z')
    ax7.set_ylabel('Precursor mass error, ppm')

    # 8 Scatter plot of ppm mass error vs retention time
    plt.subplot(grid[1, 1:3])
    ax8 = sns.lineplot(fresh_data['RT array'],fresh_data['ppm error array'],
                          alpha=0.4, s=10,c=('#A9391E',),marker='X')
    ax8.set_xlabel('RT, min')
    ax8.set_ylabel('Precursor mass error, ppm')

    # 9 Line plot of mean chromatographic peak width vs Search ID
    plt.subplot(grid[1, 3])
    ax9 = sns.lineplot(x='search_id', y='mean_peak_width',
                       data=df)
    ax9.set_xlabel('search_id')
    ax9.set_ylabel('Mean chrom peak width, s')

    # 10 Line plot of standard deviation of chromatographic peak width vs Search ID
    plt.subplot(grid[1, 4])
    ax10 = sns.lineplot(x='search_id', y='peak_width_stdev',
                       data=df)
    ax10.set_xlabel('search_id')
    ax10.set_ylabel('Chrom peak width standard dev, s')

    # 11 Histogram of the search engine PSM scores
    plt.subplot(grid[1, 5])
    # Use the hoisto with fixed number of bins = 50
    show_histo_fixedbinnum(fresh_data['Peak width array'],
                           'Chrom peak width, s',50)

    # 12 Line plot of the mean and median precursor mass error vs Search ID
    df_mzerrs = pd.melt(df,id_vars=('search_id',),
                      value_vars=('mean_mz_err_ppm', 'median_mz_err_ppm'),
                      var_name='Parameter',value_name='M/m error, ppm')
    plt.subplot(grid[2, 1])
    ax12 = sns.lineplot(x='search_id', y='M/m error, ppm',
                       hue='Parameter',data=df_mzerrs)
    ax12.set_xlabel('search_id')
    ax12.set_ylabel('Mean/median precursor mass error, ppm')
    ax12.legend(loc='lower left', fontsize='small', framealpha=0.5)

    # 13 Line plot of the mean precursor mass error standard deviation vs Search ID
    plt.subplot(grid[2, 2])
    ax13 = sns.lineplot(x='search_id', y='mz_err_ppm_stdev',
                       data=df)
    ax13.set_xlabel('search_id')
    ax13.set_ylabel('Stand deviation of mass error, ppm')

    # 14 Line plot of mean precursor intensity vs Search ID
    plt.subplot(grid[2, 3])
    ax14 = sns.lineplot(x='search_id', y='mean_prec_intensity',
                       data=df)
    ax14.set_xlabel('search_id')
    ax14.set_ylabel('Mean precursor intensity')
    
    # 15 Histogram of the Log10-transformed precursor abundances
    plt.subplot(grid[2, 4])
    # Use 0.2 to get a suitable number of bins for Log10 precursor abundance
    show_histo(np.log10(fresh_data['Prec int array']),'Log10 Precursor Intensity',0.2)  

    # 16 Line plot of the mean search engine scores vs Search ID
    plt.subplot(grid[2, 5])
    ax16 = sns.lineplot(x='search_id', y='mean_sengine_score',
                        data=df)
    ax16.set_xlabel('search_id')
    ax16.set_ylabel('Mean search engine score')
    
    # 17 Histogram of the search engine PSM scores
    plt.subplot(grid[2, 6])
    # Use the hoisto with fixed number of bins = 20
    show_histo_fixedbinnum(fresh_data['Search engine score array'],
                           'Search engine score',20)

    if onlysave is False:
        plt.show()
    elif onlysave is True:
        plt.savefig(plotsavingpath)
        print('Saved the qc plots into the file',plotsavingpath)
    # ---------------------UNCOMMENT TO SHOW THE PLOT-------------------------
        os.startfile(plotsavingpath, 'open')
    return True

def testing_load_example_files():
    df_dict = {}
    df_dict['Proteins'] = pd.read_csv('TargetProtein.txt',delimiter='\t')
    df_dict['Peptide Groups'] = pd.read_csv('TargetPeptideGroup.txt',delimiter='\t')
    df_dict['PSMs'] = pd.read_csv('TargetPeptideSpectrumMatch.txt',delimiter='\t')
    df_dict['MS/MS Spectrum Info'] = pd.read_csv('MSnSpectrumInfo.txt',delimiter='\t')
    df_dict['Input Files'] = pd.read_csv('WorkflowInputFile.txt',delimiter='\t')
    
    return df_dict

def write_new_results(conn, table, result_tuple):
    '''You need to make sure that the input tuple has the proper length and format.
        The function will return False if the writing results in error.
        But if the writing and successful, but formatting is wrong, it will not check for it.'''
    writing_successful = False
    cur = conn.cursor()

    sql_command = ("INSERT INTO " + table +
                   " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
    try:
        cur.execute(sql_command, result_tuple)
        conn.commit()
        writing_successful = True
        print('Successfuly saved the summary data into the database:')
        print(str(result_tuple))
        print()
    except:
        print('Could not write the results to the database into table',table)
        print(result_tuple)
    cur.close()
    return writing_successful

if __name__ == '__main__':

    db_files = {'Lumos':'Z:\\Lumos\\QC\\QC_Reports\\QC_DB\\qc_lumos_st191029.db',
                'Lumos_FAIMS':'Z:\\Lumos\\QC\\QC_Reports\\QC_DB\\qc_lumos_st191029.db',
                'Fusion':'Z:\\Fusion\\QC\\QC_Reports\\QC_DB\\qc_fusion_st191029.db',
                'Fusion_FAIMS':'Z:\\Fusion\\QC\\QC_Reports\\QC_DB\\qc_fusion_st191029.db',
                'QEHF':'Z:\\QExactive HF\\QC\\QC_Reports\\QC_DB\\qc_qehf_st191107.db',
                'Exp':'Z:\\Exploris\\QC\\Reports\\DB\\qc_exploris_st210914.db',
                'Exp_FAIMS':'Z:\\Exploris\\QC\\Reports\\DB\\qc_exploris_st210914.db',
                'Eclipse':'Z:\\Eclipse\\QC\\QC_Reports\\QC_DB\\qc_eclipse_st220817.db',
                'Eclipse_FAIMS':'Z:\\Eclipse\\QC\\QC_Reports\\QC_DB\\qc_eclipse_st220817.db'}
    
    index_col = 'search_id'
    colnames = ('search_id','raw_file','file_date','search_date','instrument',
                'protein_number','peptide_number','psm_number',
                'msms_number','id_rate',
                'mean_psm_it_ms','median_psm_it_ms',
                'mean_msms_it_ms','median_msms_it_ms',
                'mean_mz_err_ppm','median_mz_err_ppm','mz_err_ppm_stdev',
                'total_prec_intensity','mean_prec_intensity',
                'mean_sengine_score',
                'mean_peak_width','peak_width_stdev',
                'pept_416','pept_425','pept_488','pept_495',
                'pept_567','pept_652','pept_655',
                'comment')
    num = 15

    #-----------------------------------------@@@@-----------------------------------------@@@@
    
##    in_dfs = testing_load_example_files()

    input_json = get_console_arg()

    # Read the parameters from json files and load the PD tables into dataframes
    input_params = ius.NodeArgs(input_json)
    tableobj = ius.InputTables(input_params.return_all_table_properties())
    in_dfs = tableobj.return_all_tables()

##    for i in in_dfs.keys():
##        print(('<' + i + '>'))
##        print(in_dfs[i].columns)
        
    # The dictionary contains extracted information from the latest Proteome Discoverer search
    fresh_data = generate_input_dict(in_dfs)
    # Select the table based on the Instrument name from the search data
    table = fresh_data['Tuple for database'][4].lower()
    dbase = db_files[fresh_data['Tuple for database'][4]]
    print('Chosen the table', table)
    
    conn = sqlite3.connect(dbase)
    
##    result_tuple = generate_simulated_tuple()
##    writing_succeeded = write_new_results(conn, table, result_tuple)
    writing_succeeded = write_new_results(conn, table, fresh_data['Tuple for database'])

    df = return_last_rows(conn, table, index_col, num, colnames)

    conn.close()
    
    report_paths = create_reporting_fnames(dbase, num, fresh_data['Tuple for database'][4],
                                           fresh_data['Tuple for database'][1])
    try:
        show_plots(df,fresh_data,plotsavingpath=report_paths['graph_out'],onlysave=True)
    except:
        warnings.warn('Could not save or show the result graphs', UserWarning)
    
    try:
        df.to_excel(report_paths['excel_out'])
    except:
        warnings.warn('Could not save the Excel file', UserWarning)


    
