# Python 3.8
# write the individual topology files from the template.
# SPC_E should be SPC/E, but we can't use a backslash in a file name

models = ['SPC','SPC_E','TIP3P','TIP4P','TIP4P-EW','TIP4P_2005','OPC3','OPC','TIP4P-ST','TIP3P-ST','TIP4P-FB','GOPAL']

names = dict()
params = dict()

#open the data file, and find the parameters for each mode.
data = open("Force_field_parameters.csv", "r")
ds = data.readlines()
for d in ds:
    vals = d.split(',')
    if vals[0] == 'Lsolv_Top_variable':
        for i,v in enumerate(vals):
            if 'lam' in v:
                names[i] = v
                
    if vals[0] in models:
        parm = dict()
        for i,v in enumerate(vals):
            if i in names:
                parm[names[i]] = vals[i]
        params[vals[0]] = parm

# write out a template for each model.
for m in models:
    template = open("Lsolv_template.top", "r")
    ts = template.readlines()
    template.close()
    
    filename = m+'.top'
    f = open(filename, "w")

    for t in ts:
        for n in names:
            t = t.replace(names[n],params[m][names[n]])
        f.write(t)
    f.close()
