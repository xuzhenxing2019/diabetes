# the code is for generating diabetes_incidence
import pandas as pd
import pickle as pkl

BASELINE_LEFT = -1095 # day
BASELINE_RIGHT = 0 # day
FOLLOWUP_LEFT = 31 # day
FOLLOWUP_RIGHT = 180 # day

flag_exclusion_duration = '0_3_years' # '-7-3years','0-3years'
HbA1C_threshold = 6.5

def func_set_diabetes_flag(candidate_df):
    add_flag_baseline = ['flag_baseline_diabetes',
                         'flag_baseline_E10', 'flag_baseline_E11', 'flag_baseline_E080913','flag_baseline_O24',
                         'flag_baseline_insulin','flag_baseline_insulin_GLP_1','flag_baseline_pramlintide','flag_baseline_oral_hypoglycemic', 'flag_baseline_HbA1C']  # flag is used to exclude diabetes at baseline

    add_flag_followup = ['flag_followup_diabetes',
                         'flag_diabetes_dx', 'flag_diabetes_dx_t2e', 'flag_diabetes_med', 'flag_diabetes_med_t2e',
                         'flag_diabetes_lab', 'flag_diabetes_lab_t2e','flag_diabetes_dx_med', 'flag_diabetes_dx_med_t2e',
                         'flag_diabetes_dx_lab', 'flag_diabetes_dx_lab_t2e','flag_diabetes_med_lab', 'flag_diabetes_med_lab_t2e'
                         ]  # flag is used to include diabetes at followup based on different criteria and their time to event

    candidate_df[add_flag_baseline] = candidate_df.apply(lambda x: (0, 0, 0, 0, 0, 0, 0, 0, 0, 0), axis=1, result_type='expand')
    candidate_df[add_flag_followup] = candidate_df.apply(lambda x: (0, 0, 9999, 0, 9999, 0, 9999, 0, 9999, 0, 9999, 0, 9999), axis=1,result_type='expand')  # set 9999 as default time to event
    candidate_df['flag_diabetes'] = 0  # flag is used for final diabetes based on baseline diabetes flag and followup diabetes flag
    candidate_df['flag_diabetes_t2e'] = 9999 # time to event based on min (flag_diabetes_dx_t2e, flag_diabetes_med_t2e, flag_diabetes_lab_t2e, flag_diabetes_dx_med_t2e, flag_diabetes_dx_lab_t2e, flag_diabetes_med_lab_t2e) and maxfollowup, FOLLOWUP_RIGHT
    return candidate_df

def func_read_ori_data(candidate_folder_path, candidate_site):
    # read medication, diagnosis, lab files;
    covid_data = pd.read_csv(candidate_folder_path +candidate_site + '\\'+'matrix_cohorts_covid_for_diabetes_' + candidate_site + '.csv',index_col=False)
    covid_data['index date'] = pd.to_datetime(covid_data['index date'])

    id_medication_f = open(candidate_folder_path +candidate_site + '\\' + "medication_" + candidate_site + ".pkl", "rb")
    id_medication = pkl.load(id_medication_f)  # id_medication is dictionary, and patid is key; value is a list with each item: such as (Timestamp('2006-12-22 00:00:00'), '757969', 364)

    id_diagnosis_f = open(candidate_folder_path +candidate_site + '\\' + "diagnosis_" + candidate_site + ".pkl", "rb")
    id_diagnosis = pkl.load(id_diagnosis_f) # id_diagnosis is dictionary, and patid is key; value is a list with each item: such assample: (Timestamp('2020-08-14 00:00:00'), 'R07.2', '10', 'AV'),

    id_lab = pd.read_csv(candidate_folder_path +candidate_site + '\\' + "lab_HbA1c_" + candidate_site + ".csv", index_col=False)
    id_lab = id_lab[['patid', 'encounterid', 'lab_loinc', 'specimen_date', 'result_num', 'result_unit']]
    id_lab = id_lab.rename(columns={'patid': 'PATID', 'encounterid': 'ENCOUNTERID', 'lab_loinc': 'LAB_LOINC',
                                    'specimen_date': 'SPECIMEN_DATE', 'result_num': 'RESULT_NUM', 'result_unit': 'RESULT_UNIT'})

    return covid_data, id_medication, id_diagnosis, id_lab

def _is_in_baseline(event_time, index_time):
    return BASELINE_LEFT <= (event_time - index_time).days <= BASELINE_RIGHT

def _is_in_followup(event_time, index_time):
    return FOLLOWUP_LEFT <= (event_time - index_time).days <= FOLLOWUP_RIGHT

def func_find_HbA1C_value_date(candidate_loinc_code, candidate_loinc_value, candidate_loinc_date, HbA1C_loinc):
    value_HbA1C = []
    date_HbA1C = []
    for i in range(len(candidate_loinc_code)):
        if (candidate_loinc_code[i] in HbA1C_loinc) and (candidate_loinc_value[i] not in ['.']): # ['.'] noise data
            if (float(candidate_loinc_value[i]) >0) and (float(candidate_loinc_value[i]) <=100): # HbA1C value in [0,100] (%)
                if (float(candidate_loinc_value[i]) >= HbA1C_threshold):
                    value_HbA1C.append(float(candidate_loinc_value[i]))
                    date_HbA1C.append(candidate_loinc_date[i])

    return value_HbA1C, date_HbA1C

def func_find_specific_medication_code_date(candidate_med_code_list, candidate_med_date_list, specific_medication_rxcui_list):
    specific_rxcui = specific_medication_rxcui_list
    object_med_code = []
    object_med_date = []
    for i in range(len(candidate_med_code_list)):
        if (int(candidate_med_code_list[i]) in specific_rxcui):
            object_med_code.append(candidate_med_code_list[i])
            object_med_date.append(candidate_med_date_list[i])
    return object_med_code, object_med_date

def func_find_specific_diagnosis_code_date(candidate_dx_code_list, candidate_dx_date_list, candidate_type_list, specific_dx_list):
    Inpatient_IP_type = ['EI','IP','OS'] # inpatient
    Outpatient_OP_type = ['AV','OA','TH'] + ['ED'] # outpatient and emergency
    IP_code, OP_code = [],[]
    IP_date, OP_date = [],[]

    for i in range(len(candidate_dx_code_list)):
        # print('diabetes_code',i)
        if candidate_dx_code_list[i] in specific_dx_list:
            if candidate_type_list[i] in Inpatient_IP_type:  # for inpatients
                IP_code.append(candidate_dx_code_list[i])
                IP_date.append(candidate_dx_date_list[i])

            if candidate_type_list[i] in Outpatient_OP_type: # for outpatients
                OP_code.append(candidate_dx_code_list[i])
                OP_date.append(candidate_dx_date_list[i])

    return IP_code, IP_date, OP_code, OP_date

def _eligibility_dx_only(id_indexrecord, id_dx, func_is_in_followup, dx_CL2_list):
    print("Step: applying _eligibility_followup_dx", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]

        v_dx = id_dx.get(str(pid), []) # # checking all diagnosis codes of a patient
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx:
                dx_date, dx_code, dx_type = r[0], r[1], r[3]
                if func_is_in_followup(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)

            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)

                OP_date_list = list(set(OP_date_list)) # remove, same time, for checking E11* at different days
                OP_date_list.sort() # rank

                if IP_code_list:  # at least 1 IP E11*
                    id_indexrecord.loc[i, 'flag_diabetes_dx'] = 1
                    id_indexrecord.loc[i,'flag_diabetes_dx_t2e'] = (min(IP_date_list) - index_date).days # add time to event

                if len(OP_date_list)>=2:# OP_date_list has ranked, at least 2 OP E11* at on different days
                    id_indexrecord.loc[i, 'flag_diabetes_dx'] = 1
                    id_indexrecord.loc[i, 'flag_diabetes_dx_t2e'] = (OP_date_list[1] - index_date).days # add time to event, Note that, OP_date_list is ranked, and use the second date as time to event

                if IP_code_list and (len(set(OP_date_list))>=2):
                    id_indexrecord.loc[i, 'flag_diabetes_dx'] = 1
                    id_indexrecord.loc[i, 'flag_diabetes_dx_t2e'] = min([(min(IP_date_list) - index_date).days, (OP_date_list[1] - index_date).days])  # add min time of IP and OP as time to event
    # Summary
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_dx']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def _eligibility_med_only(id_indexrecord, id_med, id_dx, id_lab, func_is_in_followup, dx_CL2_list, lab_HbA1C_loinc_list, med_CL5f_rxcui_list, med_CL5g_rxcui_list):
    print("Step: applying _eligibility_followup_med_only",'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]

        v_dx = id_dx.get(str(pid), []) # # checking T2DM dx code: num_dx_CL2
        num_dx_CL2 = 0
        if v_dx:
            dx_date_list = []
            dx_code_list = []
            dx_type_list = []
            for r in v_dx:
                dx_date = r[0]
                dx_code = r[1]
                dx_type = r[3]
                if func_is_in_followup(dx_date, index_date):
                    dx_date_list.append(dx_date)
                    dx_code_list.append(dx_code)
                    dx_type_list.append(dx_type)
            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)
                num_dx_CL2 = len(IP_code_list)

        # checking value ---HbA1C  >=6.5 during followup period
        v_lab_df = id_lab.loc[id_lab['PATID']==pid]
        num_lab_CL6 = 0
        if not v_lab_df.empty:
            lab_loinc_code = []
            lab_loinc_value = []
            lab_loinc_date = []

            for j, row_j in v_lab_df.iterrows():
                lab_code = row_j[2]
                lab_date = pd.Timestamp(row_j[3])
                lab_value = row_j[4]
                if func_is_in_followup(lab_date, index_date):
                    lab_loinc_code.append(lab_code)
                    lab_loinc_value.append(lab_value)
                    lab_loinc_date.append(lab_date)

            if lab_loinc_code:
                HbA1C_value_list, HbA1C_date = func_find_HbA1C_value_date(lab_loinc_code, lab_loinc_value, lab_loinc_date, lab_HbA1C_loinc_list)  # find HbA1C value (>=6.5) and their date
                num_lab_CL6 = len(HbA1C_value_list)

        # considering CL5f, 5g prescription
        v_med = id_med.get(str(pid), [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                med_CL5f_list, med_CL5f_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5f_rxcui_list)
                med_CL5g_list, med_CL5g_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5g_rxcui_list)
                flag_med_CL5f = False
                flag_med_CL5g = False
                if len(med_CL5f_list) >= 1: # only CL5f criteria
                    flag_med_CL5f = True
                    id_indexrecord.loc[i, 'flag_diabetes_med'] = 1
                    id_indexrecord.loc[i, 'flag_diabetes_med_t2e'] = (min(med_CL5f_date) - index_date).days  # add time to event

                if len(med_CL5g_list) >= 1: # CL5g with T2DM or HbA1C criteria
                    sum_CL2_CL6 = num_dx_CL2 + num_lab_CL6
                    if sum_CL2_CL6 >=1:
                        flag_med_CL5g = True
                        id_indexrecord.loc[i, 'flag_diabetes_med'] = 1
                        id_indexrecord.loc[i, 'flag_diabetes_med_t2e'] = (min(med_CL5g_date) - index_date).days  # add time to event

                if flag_med_CL5f and flag_med_CL5g:
                    id_indexrecord.loc[i, 'flag_diabetes_med'] = 1
                    id_indexrecord.loc[i, 'flag_diabetes_med_t2e'] = min([(min(med_CL5f_date) - index_date).days, (min(med_CL5g_date) - index_date).days])   # add min time to event between CL5f and CL5g
    # Summary
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_med']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_med']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_med'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def _eligibility_lab_only(id_indexrecord, id_lab, id_med, func_is_in_followup, lab_HbA1C_loinc_list, med_insulin_rxcui_list):
    print("Step: applying _eligibility_followup_lab_only",'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]

        # considering insulin prescription, # set tag for flag_followup_insulin,
        flag_followup_insulin = False # if there are insulin, it is set as True
        v_med = id_med.get(str(pid), [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                med_insulin_list, med_insulin_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_insulin_rxcui_list)
                if len(med_insulin_list) > 0:
                    flag_followup_insulin = True

        # considering lab ---HbA1C
        v_lab_df = id_lab.loc[id_lab['PATID']==pid]
        if not v_lab_df.empty:
            lab_loinc_code = []
            lab_loinc_value = []
            lab_loinc_date = []

            for j, row_j in v_lab_df.iterrows():
                lab_code = row_j[2]
                lab_date = pd.Timestamp(row_j[3])
                lab_value = row_j[4]

                if func_is_in_followup(lab_date, index_date):
                    lab_loinc_code.append(lab_code)
                    lab_loinc_value.append(lab_value)
                    lab_loinc_date.append(lab_date)

            if lab_loinc_code:
                HbA1C_value_list, HbA1C_date_list = func_find_HbA1C_value_date(lab_loinc_code, lab_loinc_value, lab_loinc_date, lab_HbA1C_loinc_list) # find HbA1C value (>=6.5) and their date
                if (not flag_followup_insulin) and HbA1C_value_list: # if there are no insulin, and HbA1C >=6.5
                    id_indexrecord.loc[i, 'flag_diabetes_lab'] = 1
                    id_indexrecord.loc[i, 'flag_diabetes_lab_t2e'] = (min(HbA1C_date_list) - index_date).days  # add time to event

    # Applying excluding:
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_lab']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_lab']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_lab'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def _eligibility_dx_med(id_indexrecord, id_dx, id_med, func_is_in_followup, dx_CL2_list, med_CL5abcde_rxcui_list):
    print("Step: applying _eligibility_followup_Diagnosis_Meds", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]

        flag_followup_dx = False
        flag_followup_med = False
        # checking dx code
        v_dx = id_dx.get(str(pid), [])
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx:
                dx_date = r[0]
                dx_code = r[1]
                dx_type = r[3] # type
                if func_is_in_followup(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)
            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)
                if len(IP_code_list + OP_code_list) >= 1: # at least one E11*
                    flag_followup_dx = True
                    candidate_dx_date = min(IP_date_list+OP_date_list) #

        # checking medication
        v_med = id_med.get(str(pid), [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                med_CL5abcde_list, med_CL5abcde_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5abcde_rxcui_list)
                if len(med_CL5abcde_list) >=1:
                    flag_followup_med = True
                    candidate_med_date = min(med_CL5abcde_date)

        if flag_followup_dx and flag_followup_med:
            id_indexrecord.loc[i, 'flag_diabetes_dx_med'] = 1 # add diabetes flag
            id_indexrecord.loc[i, 'flag_diabetes_dx_med_t2e'] = (max([candidate_dx_date, candidate_med_date]) - index_date).days  # add time to event

    # Summary:
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_dx_med']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx_med']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx_med'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def _eligibility_dx_lab(id_indexrecord, id_dx, id_lab, func_is_in_followup, dx_CL2_list, lab_HbA1C_loinc_list):
    print("Step: applying _eligibility_f_Diagnosis_Lab", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]

        # set flags for followup dx, followup lab
        flag_followup_dx = False
        flag_followup_lab = False
        # checking dx code
        v_dx = id_dx.get(str(pid), [])
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx:
                dx_date = r[0]
                dx_code = r[1]
                dx_type = r[3]
                if func_is_in_followup(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)
            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)
                if len(IP_code_list + OP_code_list) >= 1: # at least one E11*
                    flag_followup_dx = True
                    candidate_dx_date = min(IP_date_list + OP_date_list)

        # checking lab ---HbA1C
        v_lab_df = id_lab.loc[id_lab['PATID']==pid]
        if not v_lab_df.empty:
            lab_loinc_code = []
            lab_loinc_value = []
            lab_loinc_date = []
            for j, row_j in v_lab_df.iterrows():
                lab_code = row_j[2]
                lab_date = pd.Timestamp(row_j[3])
                lab_value = row_j[4]
                if func_is_in_followup(lab_date, index_date):
                    lab_loinc_code.append(lab_code)
                    lab_loinc_value.append(lab_value)
                    lab_loinc_date.append(lab_date)
            if lab_loinc_code:
                HbA1C_value_list, HbA1C_date_list = func_find_HbA1C_value_date(lab_loinc_code, lab_loinc_value, lab_loinc_date, lab_HbA1C_loinc_list)  # find HbA1C value (>=6.5) and their date
                if HbA1C_value_list: # at least one HbA1C >=6.5
                    flag_followup_lab = True
                    candidate_lab_date = min(HbA1C_date_list)

        # combine dx and lab
        if flag_followup_dx and flag_followup_lab:
            id_indexrecord.loc[i, 'flag_diabetes_dx_lab'] = 1 # add diabetes flag
            id_indexrecord.loc[i, 'flag_diabetes_dx_lab_t2e'] = (max([candidate_dx_date,candidate_lab_date]) - index_date).days  # add time to event

    # Summary:
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_dx_lab']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx_lab']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_dx_lab'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def _eligibility_med_lab(id_indexrecord, id_med, id_lab, func_is_in_followup, med_CL5e_rxcui_list, med_CL5b_rxcui_list, med_CL5a_rxcui_list, lab_HbA1C_loinc_list):
    print("Step: applying _eligibility_followup_medication_lab", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, covid_flag, index_date = row[1],  row[3], row[4]
        # set three flags for flag_followup_med_lab, flag_followup_med, flag_followup_lab
        flag_followup = False
        flag_followup_med = False
        flag_followup_lab = False

        # checking medication
        v_med = id_med.get(str(pid), [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)

            if med_code_list:
                med_CL5e_list, med_CL5e_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5e_rxcui_list) #  CL 5e
                med_CL5b_list, med_CL5b_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5b_rxcui_list)  #  CL 5b
                med_CL5a_list, med_CL5a_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5a_rxcui_list)  #  CL 5a

                flag_5e = False
                flag_5b = False
                if med_CL5e_list:
                    flag_5e = True
                    flag_followup_med = True
                    candidate_med_date = min(med_CL5e_date)

                if med_CL5b_list and (not med_CL5a_list): # CL5b, without CL 5a
                    flag_5b = True
                    flag_followup_med = True
                    candidate_med_date = min(med_CL5b_date)

                if flag_5e and flag_5b:
                    candidate_med_date = min(med_CL5e_date + med_CL5b_date)

        # checking lab ---HbA1C
        v_lab_df = id_lab.loc[id_lab['PATID']==pid]
        if not v_lab_df.empty:
            lab_loinc_code = []
            lab_loinc_value = []
            lab_loinc_date = []
            for j, row_j in v_lab_df.iterrows():
                lab_code = row_j[2]
                lab_date = pd.Timestamp(row_j[3])
                lab_value = row_j[4]

                if func_is_in_followup(lab_date, index_date):
                    lab_loinc_code.append(lab_code)
                    lab_loinc_value.append(lab_value)
                    lab_loinc_date.append(lab_date)

            if lab_loinc_code:
                HbA1C_value_list, HbA1C_date_list = func_find_HbA1C_value_date(lab_loinc_code, lab_loinc_value, lab_loinc_date, lab_HbA1C_loinc_list)  # find HbA1C value (>=6.5) and their date
                if HbA1C_value_list: # at least one HbA1C >=6.5
                    flag_followup_lab = True
                    candidate_lab_date = min(HbA1C_date_list)

        if flag_followup_med and flag_followup_lab:
            id_indexrecord.loc[i, 'flag_diabetes_med_lab'] = 1 # add diabetes flag
            id_indexrecord.loc[i, 'flag_diabetes_med_lab_t2e'] = (max([candidate_med_date,candidate_lab_date]) - index_date).days  # add time to event

    # Summary:
    total_diabetes = id_indexrecord.loc[id_indexrecord['flag_diabetes_med_lab']==1]
    covid_pos_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_med_lab']==1)&(id_indexrecord['covid']==1)]
    covid_neg_diabetes = id_indexrecord.loc[(id_indexrecord['flag_diabetes_med_lab'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-diabetes: {}\tcovid pos diabetes: {}\tcovid neg diabetes: {}'.format(len(total_diabetes), len(covid_pos_diabetes), len(covid_neg_diabetes)))

    return id_indexrecord

def label_diabetes_baseline(id_indexrecord, id_dx, id_med, id_lab, func_is_in_baseline, dx_CL1_list, dx_CL2_list, dx_CL3_list, dx_CL4_list, med_CL5a_rxcui_list, med_CL5c_rxcui_list, med_CL5d_rxcui_list, med_CL5e_rxcui_list, lab_HbA1C_loinc_list):
    print("Step: checking  E10*; E11*; E08*,E09*,E13*; O24.0*,O24.1*,O24.3*,O24.8*,O24.9*, Insulin (CL 5a), Insulin/GLP-1 (CL 5c), Pramlintide (CL 5d), or oral hypoglycemic (CL 5e), HbA1c at baseline")
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid, index_date = row[1], row[4]

        ## checking the dx E10*; E08*,E09*,E13*; O24.0*,O24.1*,O24.3*,O24.8*,O24.9*
        v_dx = id_dx.get(str(pid), [])
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx:
                dx_date = r[0]
                dx_code = r[1]
                dx_type = r[3]
                if func_is_in_baseline(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)
            if dx_code_list:
                E10_IP_code_list, E10_IP_date_list, E10_OP_code_list, E10_OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL1_list)
                E11_IP_code_list, E11_IP_date_list, E11_OP_code_list, E11_OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)
                E080913_IP_code_list, E080913_IP_date_list, E080913_OP_code_list, E080913_OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL3_list)
                EO24_IP_code_list, EO24_IP_date_list, EO24_OP_code_list, EO24_OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL4_list)

                id_indexrecord.loc[i, 'flag_baseline_E10'] = min([len(E10_IP_code_list + E10_OP_code_list),1])
                id_indexrecord.loc[i, 'flag_baseline_E11'] = min([len(E11_IP_code_list + E11_OP_code_list),1])
                id_indexrecord.loc[i, 'flag_baseline_E080913'] = min([len(E080913_IP_code_list + E080913_OP_code_list),1])
                id_indexrecord.loc[i, 'flag_baseline_O24'] = min([len(EO24_IP_code_list + EO24_OP_code_list),1])

        ## checking medication (Insulin (CL 5a), Insulin/GLP-1 (CL 5c), Pramlintide (CL 5d), or oral hypoglycemic (CL 5e)  )
        v_med = id_med.get(str(pid), [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_baseline(med_date, index_date):
                    med_date_list.append(med_date)
                    med_code_list.append(med_code)
            if med_code_list:
                med_CL5a_insulin, med_CL5a_insulin_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5a_rxcui_list)
                med_CL5c_insulin_GLP_1, med_CL5c_insulin_GLP_1_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5c_rxcui_list)
                med_CL5d_pramlintide, med_CL5d_pramlintide_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5d_rxcui_list)
                med_CL5e_oral_hypoglycemic, med_CL5e_oral_hypoglycemic_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_CL5e_rxcui_list)

                id_indexrecord.loc[i,'flag_baseline_insulin'] = min([len(med_CL5a_insulin),1])
                id_indexrecord.loc[i, 'flag_baseline_insulin_GLP_1'] = min([len(med_CL5c_insulin_GLP_1),1])
                id_indexrecord.loc[i, 'flag_baseline_pramlintide'] = min([len(med_CL5d_pramlintide),1])
                id_indexrecord.loc[i, 'flag_baseline_oral_hypoglycemic'] = min([len(med_CL5e_oral_hypoglycemic),1])

        # checking the value of HbA1C >=6.5 at baseline
        v_lab_df = id_lab.loc[id_lab['PATID']==pid]
        if not v_lab_df.empty:
            lab_loinc_code = []
            lab_loinc_date = []
            lab_loinc_value = []
            for j, row_j in v_lab_df.iterrows():
                lab_code = row_j[2]
                lab_date = pd.Timestamp(row_j[3])
                lab_value = row_j[4]
                if func_is_in_baseline(lab_date, index_date):
                    lab_loinc_code.append(lab_code)
                    lab_loinc_date.append(lab_date)
                    lab_loinc_value.append(lab_value)
            if lab_loinc_code:
                HbA1C_value_list, HbA1C_date_list = func_find_HbA1C_value_date(lab_loinc_code, lab_loinc_value, lab_loinc_date, lab_HbA1C_loinc_list)  # find HbA1C value (>=6.5) and their date
                id_indexrecord.loc[i, 'flag_baseline_HbA1C'] = min([len(HbA1C_value_list),1])

    return id_indexrecord

def determine_time2event(id_indexrecord):
    print("Step: adding time to event:", len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # Unnamed: 0,patid,site,covid,index date,hospitalized,ventilation,criticalcare,maxfollowup,death
        pid = row[1]
        maxfollowup = id_indexrecord.loc[i,'maxfollowup']
        flag_diabetes = id_indexrecord.loc[i,'flag_diabetes']
        flag_diabetes_dx_t2e = id_indexrecord.loc[i,'flag_diabetes_dx_t2e']
        flag_diabetes_med_t2e = id_indexrecord.loc[i,'flag_diabetes_med_t2e']
        flag_diabetes_lab_t2e = id_indexrecord.loc[i,'flag_diabetes_lab_t2e']
        flag_diabetes_dx_med_t2e = id_indexrecord.loc[i,'flag_diabetes_dx_med_t2e']
        flag_diabetes_dx_lab_t2e = id_indexrecord.loc[i,'flag_diabetes_dx_lab_t2e']
        flag_diabetes_med_lab_t2e = id_indexrecord.loc[i,'flag_diabetes_med_lab_t2e']

        if flag_diabetes ==1:
            id_indexrecord.loc[i, 'flag_diabetes_t2e'] = min([flag_diabetes_dx_t2e, flag_diabetes_med_t2e, flag_diabetes_lab_t2e, flag_diabetes_dx_med_t2e, flag_diabetes_dx_lab_t2e, flag_diabetes_med_lab_t2e])
        else:
            id_indexrecord.loc[i, 'flag_diabetes_t2e'] = min([maxfollowup, FOLLOWUP_RIGHT])

    return id_indexrecord


if __name__ == '__main__':
    COVID_data_path = 'D:\zhenxing_xu\code\pasc\pasc_phenotype\data\V15_COVID19\input_diabetes\\'
    result_path = 'D:\Zhenxing_Xu\diabetes_incidence\\'
    diabetes_CL_file = 'D:\Zhenxing_Xu\Yongkang_pro\diabetes\RECOVER Diabetes New Incidence CP Code Lists_4.22.22FINAL_Updated 1pm.xlsx'

    diabetes_CL1 = pd.read_excel(diabetes_CL_file, sheet_name='1. E10 - Type 1 DM Dx')
    diabetes_CL2 = pd.read_excel(diabetes_CL_file, sheet_name='2. E11 - Type 2 DM Dx')
    diabetes_CL3 = pd.read_excel(diabetes_CL_file, sheet_name='3. E08|E09|E13 - Other DM Dx')
    diabetes_CL4b = pd.read_excel(diabetes_CL_file, sheet_name='4b. O24 No O24.4  - Gest DM Dx')
    diabetes_CL5a = pd.read_excel(diabetes_CL_file, sheet_name='5a. Insulin')
    diabetes_CL5b = pd.read_excel(diabetes_CL_file, sheet_name='5b. GLP-1')
    diabetes_CL5c = pd.read_excel(diabetes_CL_file, sheet_name='5c. Insulin|GLP-1')
    diabetes_CL5d = pd.read_excel(diabetes_CL_file, sheet_name='5d. Pramlintide')
    diabetes_CL5e = pd.read_excel(diabetes_CL_file, sheet_name='5e. All Oral Hypoglycemic')
    diabetes_CL5f = pd.read_excel(diabetes_CL_file, sheet_name='5f. Oral Hypo No Met|Phen|SGLT2')
    diabetes_CL5g = pd.read_excel(diabetes_CL_file, sheet_name='5g. Oral Hy Met|Phen|SGLT2 Only')
    diabetes_CL6 = pd.read_excel(diabetes_CL_file, sheet_name='6. HbA1C LOINC')

    dx_CL1 = list(diabetes_CL1['code1'])
    dx_CL2 = list(diabetes_CL2['code1'])
    dx_CL3 = list(diabetes_CL3['code1'])
    dx_CL4 = list(diabetes_CL4b['code1'])
    lab_HbA1C_loinc = list(diabetes_CL6['LOINC'])
    med_insulin_rxcui = list(diabetes_CL5a['rxcui'])
    med_GLP_rxcui = list(diabetes_CL5b['rxcui'])
    med_insulin_GLP_rxcui = list(diabetes_CL5c['rxcui'])
    med_pramlintide_rxcui = list(diabetes_CL5d['rxcui'])
    med_hypoglycemic_rxcui = list(diabetes_CL5e['rxcui'])
    med_CL5f_rxcui = list(diabetes_CL5f['rxcui'])
    med_CL5g_rxcui = list(diabetes_CL5g['rxcui'])

    med_CL5abcde_rxcui = med_insulin_rxcui + med_GLP_rxcui + med_insulin_GLP_rxcui + med_pramlintide_rxcui + med_hypoglycemic_rxcui

    # build a dataframe for saving results
    diabetes_incident_cohort = pd.DataFrame()
    site_s = ['COL'] # , 'COL', 'MSHS', 'MONTE', 'NYU', 'WCM'

    for site in site_s:
        # read covid data
        id_indexrecord, id_medication, id_diagnosis, id_lab = func_read_ori_data(COVID_data_path, site)
        id_indexrecord = func_set_diabetes_flag(id_indexrecord)
        print('site and number of patients:',site, len(id_indexrecord))

        # Type 2--diabetes--6 criteria
        ## EC--dx_only
        id_indexrecord = _eligibility_dx_only(id_indexrecord, id_diagnosis, _is_in_followup, dx_CL2)

        ## EC--med_only
        id_indexrecord = _eligibility_med_only(id_indexrecord, id_medication, id_diagnosis, id_lab, _is_in_followup, dx_CL2, lab_HbA1C_loinc, med_CL5f_rxcui, med_CL5g_rxcui)

        ## EC--lab_only
        id_indexrecord = _eligibility_lab_only(id_indexrecord, id_lab, id_medication, _is_in_followup, lab_HbA1C_loinc, med_insulin_rxcui)

        ## EC--dx--meds
        id_indexrecord = _eligibility_dx_med(id_indexrecord, id_diagnosis, id_medication, _is_in_followup, dx_CL2, med_CL5abcde_rxcui)

        ## EC--dx--lab
        id_indexrecord = _eligibility_dx_lab(id_indexrecord, id_diagnosis, id_lab, _is_in_followup, dx_CL2, lab_HbA1C_loinc)

        ## EC--med--lab
        id_indexrecord = _eligibility_med_lab(id_indexrecord, id_medication, id_lab, _is_in_followup, med_hypoglycemic_rxcui, med_GLP_rxcui, med_insulin_rxcui,lab_HbA1C_loinc)

        ## baseline exclusion
        id_indexrecord = label_diabetes_baseline(id_indexrecord, id_diagnosis, id_medication, id_lab, _is_in_baseline, dx_CL1, dx_CL2, dx_CL3, dx_CL4, med_insulin_rxcui, med_insulin_GLP_rxcui, med_pramlintide_rxcui, med_hypoglycemic_rxcui, lab_HbA1C_loinc)

        ## baseline diabetes
        id_indexrecord['flag_baseline_diabetes'] = id_indexrecord['flag_baseline_E10'] + id_indexrecord['flag_baseline_E11'] + id_indexrecord['flag_baseline_E080913'] + id_indexrecord['flag_baseline_O24'] \
                                                   + id_indexrecord['flag_baseline_insulin'] + id_indexrecord['flag_baseline_insulin_GLP_1'] + id_indexrecord['flag_baseline_pramlintide']+ id_indexrecord['flag_baseline_oral_hypoglycemic'] + id_indexrecord['flag_baseline_HbA1C']
        id_indexrecord.loc[id_indexrecord['flag_baseline_diabetes']>=1,'flag_baseline_diabetes'] = 1

        ## followup diabetes
        id_indexrecord['flag_followup_diabetes'] = id_indexrecord['flag_diabetes_dx'] + id_indexrecord['flag_diabetes_med'] + id_indexrecord['flag_diabetes_lab'] \
                                                   + id_indexrecord['flag_diabetes_dx_med'] + id_indexrecord['flag_diabetes_dx_lab'] + id_indexrecord['flag_diabetes_med_lab']
        id_indexrecord.loc[id_indexrecord['flag_followup_diabetes'] >= 1, 'flag_followup_diabetes'] = 1

        ## final diabetes based on baseline and followup diabetes
        id_indexrecord.loc[(id_indexrecord['flag_followup_diabetes']==1)&(id_indexrecord['flag_baseline_diabetes']==0),'flag_diabetes'] = 1

        ## determine time to event
        id_indexrecord = determine_time2event(id_indexrecord)

        diabetes_incident_cohort = diabetes_incident_cohort.append(id_indexrecord, ignore_index=True)
    # save final results to csv
    diabetes_incident_cohort.to_csv(result_path + 'diabetes_incident' + '_exclusion_duration_'+flag_exclusion_duration+'_0816'+'.csv', index=False)
    print('len(diabetes_incident_cohort)', len(diabetes_incident_cohort)) # 361401











