import numpy as np
import scipy.stats as stats
import csv as csv

def create_table(group,var):
    
    group = group[var!=999]
    var = var[var!=999]

    if np.min(group) > 0:
        group = group-1
    if np.min(var) > 0:
        var = var-1

    table = np.zeros((np.max(group)+1,np.max(var)+1))
    
    for i in range(0,max(group)+1):
        for j in range(0,max(var)+1):
            table[i,j] = np.sum(var[group==i]==j)
  
    table = np.flipud(table)
    
    if table.shape[1] == 2:
        table = np.fliplr(table)
    return table

def odds_ratio_calc(table):
    odds_ratio = np.divide(table[0,0],table[0,1]) / np.divide(table[1,0],table[1,1])
    x2, pvalue, _ ,_ = stats.chi2_contingency(table,correction=True)

    std_error = np.sqrt((np.divide(1,[table[0,0]]) 
                        + np.divide(1,[table[1,0]]) 
                        + np.divide(1,[table[0,1]])  
                        + np.divide(1,[table[1,1]]))
                       )
    lci = np.log(odds_ratio) - (1.96*std_error)
    uci = np.log(odds_ratio) + (1.96*std_error)

    lci = np.exp(lci)
    uci = np.exp(uci)
    return odds_ratio, pvalue, lci[0], uci[0], x2

def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.nanmean(x) - np.nanmean(y)) / np.sqrt(((nx-1)*np.nanstd(x, ddof=1) ** 2 + (ny-1)*np.nanstd(y, ddof=1) ** 2) / dof)

def ptsd_ttest(filename,data,var_labels,PTSD):
    
    for i in range(0,data.shape[1]):
        group = PTSD
        _data = data[:,i]
        group = group[_data!=-1]
        _data = _data[_data!=-1]

        t_val,p_val = stats.ttest_ind(_data[group==0],_data[group==1], nan_policy='omit')

        with open(filename,'a') as newFile:
            newFileWriter = csv.writer(newFile)
            newFileWriter.writerow([var_labels[i],
                                    np.round(np.nanmean(_data[group==0]),2),
                                    np.round(np.nanstd(_data[group==0]),2),
                                    np.round(np.nanmean(_data[group==1]),2),
                                    np.round(np.nanstd(_data[group==1]),2),
                                    np.round(np.nanmean(_data),2),
                                    np.round(np.nanstd(_data),2),
                                    np.round(t_val,2),np.round(p_val,2),
                                    np.round(cohen_d(_data[group==0],_data[group==1]),2)
                                   ])
    
def print_OR(filename,table,lab,OR,uci,lci,p,x2):

    _p1 = np.multiply(np.divide(table[1,0],np.sum(table[1,:])),100)
    _p2 = np.multiply(np.divide(table[0,0],np.sum(table[0,:])),100)
    _p3 = np.sum(table[:,0])
    _p4 = np.multiply(_p3 / np.sum(table),100)

    output = [lab, table[1,0],np.round(_p1,2),table[0,0],np.round(_p2,2),
              _p3,np.round(_p4,2),np.round(OR,2), np.round(lci,2), 
              np.round(uci,2),np.round(x2,2), np.round(p,3)]
    
    with open(filename,'a') as newFile:
        newFileWriter = csv.writer(newFile)
        newFileWriter.writerow(output)
        
    if table.shape[1] > 2:
        for i in range(1,table.shape[1]):
            _p1 = np.multiply(np.divide(table[1,i],np.sum(table[1,:])),100)
            _p2 = np.multiply(np.divide(table[0,i],np.sum(table[0,:])),100)
            _p3 = np.sum(table[:,i])
            _p4 = np.multiply(_p3 / np.sum(table),100)
            
            output = [[lab,i], table[1,i],np.round(_p1,2),table[0,i],
                      np.round(_p2,2),_p3,np.round(_p4,2),np.round(OR,2), 
                      np.round(lci,2), np.round(uci,2),np.round(x2,2), 
                      np.round(p,3)
                      ]
                
            with open(filename,'a') as newFile:
                newFileWriter = csv.writer(newFile)
                newFileWriter.writerow(output)

def results_master(filename,data,label,group):
    for i in range(0,data.shape[1]):  # loop through each variable
        table = create_table(group,data[:,i]) #creates contingency table
        OR, p, lci, uci, x2 = odds_ratio_calc(table) # calculates stats.
        print_OR(filename,table,label[i],OR,uci,lci,p,x2)