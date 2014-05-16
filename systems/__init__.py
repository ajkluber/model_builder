import system

'''
Fri Mar 21 16:12:34 CDT 2014
Alexander Kluber

    This is the start of the new format for the System class. Hopefully
the cleaning style will enable it to be more.
'''

def get_system(opts):
    System = system.System(opts["Subdir"],Tf_it=opts["Tf_Iteration"],
                                Tf_act_dir=opts["Tf_Active_Directory"],Tf_refine=opts["Tf_Refinements"],
                                Mut_it=opts["Mut_Iteration"],Mut_act_dir=opts["Mut_Active_Directory"])
    return System

def new_system(subdir):
    System = system.System(subdir)
    return System

def load_systems(subdirs):
    ''' Create systems from saved options in system.info'''
    Systems = []
    for subdir in subdirs:
        print "Loading system from subdirectory: ", subdir
        System = load_system(subdir)
        Systems.append(System)
    return Systems

def new_systems(subdirs):
    ''' Create systems from saved options in system.info'''
    Systems = []
    for subdir in subdirs:
        print "Creating system for: ", subdir
        System = new_system(subdir)
        Systems.append(System)
    print "Proceeding...\n"
    return Systems

def load_system(subdir):
    ''' Given subdir that contains model.info options file. Read in options and
        create corresponding model.'''
    info_file = open(subdir+'/system.info','r')
    line = info_file.readline()
    options = {}
    while line != '':
        field = line.split()[1]
        value = info_file.readline()
        if field == "Tf_refinements":
            temp = [int(value[:-1].split()[0])]
            line = info_file.readline()
            while (line[:1] != '[') and (line != ''):
                temp.append(int(line[:-1].split()[0]))
                line = info_file.readline()
                print line
            options[field] = temp
        else:
            options[field] = value[:-1]
            line = info_file.readline()
    options = check_fields(options)
    System = get_system(options)
    print "Proceeding...\n"
    return System

def check_fields(fields):
    ''' Update names of inputted fields for backward compatibility.'''

    print "Checking system info:"
    for key in fields.keys():
        print "  ", key , " = ", fields[key]

    ## For backwards compatibility. Convert old names to new names.
    old_fields = { "Main_Path":"Path","subdir":"Subdir",
                   "Tf_iteration":"Tf_Iteration",
                   "Tf_active_directory":"Tf_Active_Directory",
                   "Tf_refinements":"Tf_Refinements",
                   "mutation_iteration":"Mut_Iteration",
                   "mutation_active_directory":"Mut_Active_Directory"}
    
    ## Replace old style keys with new style keys. 
    for key in fields.keys():
        if key in old_fields.keys():
            fields[old_fields[key]] = fields[key]
        
    fields["Tf_Iteration"] = int(fields["Tf_Iteration"])
    fields["Mut_Iteration"] = int(fields["Mut_Iteration"])
    if type(fields["Tf_Refinements"]) == str:
        fields["Tf_Refinements"] = [ int(fields["Tf_Refinements"]) ]
    elif type(fields["Tf_Refinements"]) == list:
        pass

    print "Using system info:"
    for key in fields.keys():
        print "  ", key , " = ", fields[key]

    return fields
