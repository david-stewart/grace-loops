#!/gpfs/loomis/apps/avx/software/Python/3.6.4-foss-2018a/bin/python3
from string import Template
import subprocess, os, argparse
from sys import argv
from glob import glob
from pathlib import Path
import re

''' To do:
    Add iterator for each TClonesArray object of given <name> and <class>
 [x]      int  cnt_<name>;    // how many there are in the TClonesArray for given event
 [x]      int  index_<name>;  // current index in TClonesArray
 [x]      <class>* <name>;    // pointer to index_<name> instance
 [x]      bool next_<name>(); // loops through index_<name>, setting <class>* name as it goes
 [x]      bool reset_<name>(); // resets the iterator for the class
 [x]      int  get_cnt_<name>(); // get how many there are

 [ ]      Add a class to events to match and read is###### triggers.
'''

def tag_value(f_name, tag):
    text = Path(f_name).read_text()
    parts = re.search(f'(.*//\s*:TAG START:\s*{tag}\s*)(.*)(\n\s*//\s*:TAG END:\s*{tag}\s*.*)',text,re.DOTALL)
    if not parts:
        print(f'fatal: couldn\'t find tag: ":TAG START: {tag} :TAG END: {tag}" in {f_name}')
        exit(2)
    else:
        return [ parts.group(1), parts.group(2), parts.group(3) ]

def swap_tag_value(dir_old, dir_new, f_old, tag, f_new=''):
    if not f_new:
        f_new = f_old
    '''swap the keyed section in dir_old/f_name with the keyed section in dir_new/f_name'''
    text_old = tag_value(f'{dir_old}/{f_old}',tag)
    text_new = tag_value(f'{dir_new}/{f_new}',tag)
    with open(f'{dir_old}/{f_old}','w') as fout:
        fout.write(text_old[0])
        fout.write(text_new[1])
        fout.write(text_old[2])

def set_tag_value(f_name, tag, text):
    text_old = tag_value(f_name, tag)
    with open(f_name,'w') as fout:
        fout.write(text_old[0])
        fout.write(text)
        fout.write(text_old[2])


def swap_trees(dir_new, dir_old='.'):
    swap_tag_value(f'{dir_old}/src', f'{dir_new}/src', 'events.h','Array Sizes')
    swap_tag_value(f'{dir_old}/src', f'{dir_new}/src', 'events.h','Leaf Types')
    swap_tag_value(f'{dir_old}/src', f'{dir_new}/src', 'events.h','Declare Branches')
    swap_tag_value(f'{dir_old}/src', f'{dir_new}/src', 'events.cxx','Set Branches')

def get_members(class_name, text):
    ''' 
        Allready know that name is in file_name.
        Find all the \w+; and ((\w+,)+\s*)\w+;
        inbetween class\s+<name> ... ClassDef\s*\(name,\d+\);
    '''
    # strip // comments out of the text
    strip_text = []
    for line in text.split('\n'):
        strip_text.append(line.split('//')[0])
    text = '\n'.join(strip_text)

    # parts = re.search(f'.*class\s+{class_name}(.*)ClassDef\s*\({class_name}\s*,\s*\d+\s*\)\s*;',text)
    parts = re.search(f'.*(class\s+{class_name}.*ClassDef\s*\(\s*{class_name}\s*,\s*\d+\s*\)\s*;)',text,re.DOTALL)
    if not parts:
        parts = re.search(f'.*(struct\s+{class_name}.*ClassDef\s*\(\s*{class_name}\s*,\s*\d+\s*\)\s*;)',text,re.DOTALL)
    if not parts:
        print("fatal: failed to find ClassDef for: %s in get_members"%class_name)
        exit(2)
    local_text = parts.group(1)
    # print('-----------START ',class_name)
    # print(f"Finding {class_name}'s members in\n: {local_text}")
    # print('-----------END ',class_name)

    members = []
    for x in re.finditer('\s+(\w+)\s*;',local_text):
        if not x.group(0).strip()[:2] == '//':
            # print('found member: ', x.group(1))
            members.append(x.group(1))

    # find those that end in ;        
    for x in re.finditer('\s+((\w+)\s*,\s*)+\w+\s*;',local_text):
        # print(x.group(0))
        for item in re.finditer('(\w+)\s*,',x.group(0)):
            members.append(item.group(1))
            # print('doubles: found: ',item.group(0),' |',item.groups(),'| ')
    if 'const' in members:
            members.remove('const')
    return members;


def re_leaf_wTObjects(path, *headers):
    '''Will look for src/events.h and src/events.cxx locally, and 
    in-lists/list_0.list, to find a tree of objects. It will then
    look through teh header files to see if it can find the missing classes.
    If it can, it will modify the tree not use "fChain->SetMakeClass(1);", but
    will rather generate the classes appropriately using the header files
    provided.'''

    # first, get all the output for a header file
    for file in ('src/events.h','src/events.cxx','in-lists/list_0.list'):
        if not os.path.isfile(f'{path}/{file}'):
            print(f'Fatal error, cannot find {path}/src/{file} in re_leaf_wTObjects')
            exit(2)

    h_files = []
    for file in headers:
        if "${HOME}" in file:
            file = file.replace("${HOME}",os.getenv("HOME"))
        h_files.append(file)
    
    for file in h_files:
        if not os.path.isfile(f'{file}'):
            print(f'Fatal error, cannot find {file} in re_leaf_wTObjects')
            exit(2)
    
    # TClonesArray_branches = []
    # instance_branches = []
    # regular_branches = []

    # find the missing classes
    missing_classes = dict()
    mclasses_members = dict() # all matches between class <name> and ClassDef(<name>), 
                              # that are \w+; or (\w+,)*\s+;
    in_file = open(f'{path}/in-lists/list_0.list').readline().strip()

    _out = subprocess.Popen([f'''
root -l -b {in_file}<<EOF
    // events->MakeClass("__temp")
    events->Print()
EOF
'''],shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
     )


    stdout,stderr = _out.communicate()
    try:
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
    except:
        pass

    for line in stdout.split('\n'):
        parts = re.search('Warning in .* no dictionary for class (\w+) is available',line)
        if parts:
            missing_classes[parts.group(1)] = True # switch to false when they are found
   
    
    # read in the TClonesArray values
    # for line in stdout.split('\n'):
        # print('line:',line)
    # see if the missing classes are present in the header files
    n_missing = len(missing_classes)
    req_headers = set()
    # print(h_files, n_missing)
    for h_file in h_files:
        print(h_file)
        text = Path(h_file).read_text()
        for K,V in missing_classes.items():
            if not V:
                continue
            if re.search('ClassDef\s*\\(\w*'+K,text):
                missing_classes[K] = False
                n_missing -= 1
                req_headers.add(h_file.split('/')[-1])
                # print('calling',K,'from',h_file)
                mclasses_members[K] = get_members(K,text)
                # print(mclasses_members)
                if n_missing == 0:
                    break
            else:
                print(f'NOT found {K}')
        
    # print(missing_classes)
    # print(req_headers)

    if n_missing:
        print("fatal error: still missing following classes definitions (not found in headers:")
        for (K,V) in missing_classes.items():
            if V:
                print('   ', K)
        exit(2)
    
    # OK: we know which classes are missing, and we have the required headers
    # now figure out which branches are TClonesArray of headers, 
    # which are :
    #   Regular Branches
    #   TClonesArray of type class (dict -> Members)
    #   Member of type class       (dict -> Members)
    # reg_branches = []
    # clones_branches = dict()
    # class_branches = dict()
    text = tag_value(f"{path}/src/events.h","Declare Branches")[1]
    branch_names = [ x.group(1) for x in re.finditer('\n\s+TBranch\s+\\*(\w+);',text)]
    # for x in branch_names:
        # print(x)

    leaf_types = []
    text = tag_value(f"{path}/src/events.h","Leaf Types")[1]
    for x in re.finditer('\n\s*(\S+)\s+(\S+);',text):
        if not x.group(0).strip()[:2] == '//':
            leaf_types.append((x.group(1),x.group(2)))
    
    set_branches = []
    text = tag_value(f"{path}/src/events.cxx","Set Branches")[1]
    for x in re.finditer('\n\s*fChain->SetBranchAddress\("([^"]+)",\s*([^,]*),\s*([^\)]*)',text):
        if not x.group(0).strip()[:2] == '//':
            set_branches.append((x.group(1),x.group(2),x.group(3)))

    # for L in set_branches:
        # print(L)

    # print(len(leaf_types), len(branch_names), len(set_branches))

    if not len(leaf_types) == len(branch_names):
        print('msg.01 fatal error: len(leaf_types) != len(branch_names)')
        print(' leaf types:')
        for x in leaf_types:
            print(x)
        print(' branch names:')
        for x in branch_names:
            print(x)
        exit(2)

    if not len(leaf_types) == len(set_branches):
        print('msg.02 fatal error: len(leaf_types) != len(set_branches)')
        print(' leaf types:')
        for x in leaf_types:
            print(x)
        print(' branch names:')
        for x in set_branches:
            print(x)
        exit(2)

    # ok, make a matrix that determines what each branch is
    # types : [i] -> 'normal'
    #             -> 'skip' (don't need the 'Fjet_' members, fUniqueID, or fBits)         
    #             -> ['clones', (member0, member1, member2,...)]
    #             -> ['class' , (member0, member1, member2,...)]          
    keep = len(set_branches)*['',] # keep, name (class or clones), or skip (skip clones, fBits, and fUniqueID
    current_name = ''
    current_dict = '' # clones, class, or ''
    clones  = dict()
    classes = dict()
    members = []
    i_start = 0
    fUniqueID = []; # names of fUnique branches to skip
    fBits   = []; # names of fBits branches to skip

    for i, tup in enumerate(set_branches):
        # print ('this:',tup[0][-9:],'  ',tup)
        if tup[0][-9:] == 'fUniqueID':
            fUniqueID.append(tup[0])
            # start of a class
            members = []
            if current_dict:
                current_dict[current_name] = [tuple(members),]
            if '.' in tup[0]:
                current_dict = clones
                current_name = tup[0][:-10]
                # current_type = 'clones'
                keep[i-1] = 'skip'
            else:
                current_dict = classes
                current_name = tup[2][3:-10]
            keep[i] = current_name
            # print('current_name:', current_name)
        elif current_dict is clones:
            if tup[0][:len(current_name)+1] == current_name+'.':
                keep[i] = 'skip'
                if tup[0][-5:] == 'fBits':
                    fBits.append(tup[0])
                else:
                    members.append(tup[0][len(current_name)+1:])
                    # print ("Add member",members[-1])
            else:        
                keep[i] = 'keep'
                current_dict[current_name] = [tuple(members),]
                # types[i_start].append(tuple(members))
                # members = []
                current_name = ''
                current_dict = ''
        elif current_dict is classes:
            if tup[2][3:len(current_name)+3] == current_name:
                keep[i] = 'skip'
                if tup[0][-5:] == 'fBits':
                    fBits.append(tup[0])
                else:
                    members.append(tup[0])
                    # print ("Add member",members[-1])
            else:        
                keep[i] = 'keep'
                current_dict[current_name] = [tuple(members),]
                # types[i_start].append(tuple(members))
                members = []
                current_name = ''
                current_dict = ''
        else:        
            keep[i] = 'keep'
    # get last clone or class if at the end of the tree:
    if current_dict:
        current_dict[current_name] = [tuple(members),]


    # for i in range(len(leaf_types)):
    #     print(f'{keep[i]} : {set_branches[i]}')
    #     if keep[i] in classes:
    #         print(f'class: {keep[i]} {classes[keep[i]]}')
    #     if keep[i] in clones:
    #         print(f'TClonesArray: {keep[i]} {clones[keep[i]]}')
    
    # now need to figure out which object each class is.
    # This ins't full-proof, but search for members in the class, from 
    # class <name> .... ClassDef(<name>,#);
    # If there are multiple fits, (such as Jet's members are a subset of JetwArea,
    # if both were present, then just warn the user and insert both)
    # print('-----')
    best_matches_clones  = dict()
    best_matches_classes = dict()
    for array, best_matches in zip((clones,classes),(best_matches_clones,best_matches_classes)):
        for key, members in array.items():
            matches = dict() # each class_name -> (n_missing, n_extra)
            for class_name, class_members in mclasses_members.items():
                n_missing = len(members[0])
                n_extra   = len(class_members)
                for M in members[0]:
                    if M in class_members:
                        n_missing -= 1
                        n_extra   -= 1
                matches[class_name] = (n_missing, n_extra)
            best_matches[key] = sorted(matches.items(), key=lambda x: x[1][0]+x[1][1])    
    
    # ok, ready to modify the files in place
    
    # Do Array Sizes
    text = '/* No need to use fixed array sizes because TClonesArray* are used instead of Arrays */'
    set_tag_value('src/events.h','Array Sizes',text)

    # Do "Leaf Types"
    tup = tag_value('src/events.h','Leaf Types')
    f_out = open('src/events.h','w')
    # write in the missing headers
    includes = tup[0].split('#include')
    f_out.write('#include'.join(includes[:-1]))
    f_out.write('#include')
    for header in headers:
        if header not in tup[0]:
            f_out.write(f' "{header.split("/")[-1]}"\n#include')
    f_out.write(includes[-1])

    f_out.write("//\n//  note:  TOjects have been substituted in from header files:\n")
    for header in headers:
        f_out.write(f"// {header} \n")
    f_out.write( "// The following TObjects are present in the TTree and should be used.\n")
    f_out.write( "//\n// TObjects used in TTree + members found in above headers:\n")
    # print(mclasses_members)
    for name, members in mclasses_members.items():
        f_out.write(f'//  - {name} : {" ".join(members)}\n')
    if clones: 
        f_out.write('//\n// TClonesArrays branches found in TTree (with their associated members):\n')
    for name, members in clones.items():
        f_out.write(f'// Note: a \'pseudo-iterator is provided for each TClonesArray. See leafs below.')
        f_out.write(f'//  - {name} : {" ".join(members[0])}\n')
    f_out.write('//\n// Branches of TObjects found in TTree (with their associated members):\n')
    for name, members in classes.items():
        f_out.write(f'//  - {name} : {" ".join(members[0])}\n')
    f_out.write('//\n// Best matches (according to members) of Tree branches to TObject types"\n')
    f_out.write('// Note: if any of these "first place guesses" are wrong, fix them by hand!\n')
    f_out.write('// "Branch name" -> (class-name):(#members in TTree not in header, #members in header not in TTree) ...\n')
    f_out.write('// (TCloneArray Branches)\n')
    for name, matches in best_matches_clones.items():
        f_out.write(f'//  - {name} ->')
        for entry in matches:
            f_out.write(f' {entry[0]}:{entry[1]}')
        f_out.write('\n')
    f_out.write('// (TObject Branches)\n')
    for name, matches in best_matches_classes.items():
        f_out.write(f'//  - {name} ->')
        for entry in matches:
            f_out.write(f' {entry[0]}:{entry[1]}')
        f_out.write('\n')
    # f_out.write('\n')  
    # Write the Leaf types
    for data, status in zip(leaf_types, keep):
        if status == 'keep':
            if '*' in data[1]:
                f_out.write(f'    {data[0]:19s} {data[1]} {{nullptr}};\n')
            else:
                f_out.write(f'    {data[0]:19s} {data[1]};\n')
        elif status == 'skip':
            continue
        else:
            if status in clones:
                _class = best_matches_clones[status][0][0]
                f_out.write(f'\n    //--accessor and iterator for TClonesArray* (tca) {status}\n')
                f_out.write(f'    %-19s *tca_{status} {{new TClonesArray("{_class}")}};\n'%'TClonesArray')
                f_out.write(f'    %-19s *get_{status}(int=-1);\n'%_class);
                f_out.write(f'    %-19s  iter_{status} {{ tca_{status} }};\n'%('iterTCA<%s>'%_class));
                f_out.write(f'    %-19s  size_{status}() {{ return tca_{status}->GetEntriesFast(); }};\n\n'%'int');
                # f_out.write(f'    %-19s _size_{status}{{ 0}};\n'%'int');
                # f_out.write(f'    %-19s index_{status};\n'%'int');
                # f_out.write(f'    %-19s *{status} {{nullptr}};\n'%_class);
                # f_out.write(f'    %-19s next_{status}();\n'%'bool');
                # f_out.write(f'    %-19s reset_{status}();\n'%'void');
                # f_out.write(f'    //--end {status} psuedo-iterator\n')
            elif status in classes:
                f_out.write(f'    {best_matches_classes[status][0][0]:19s} *{status};\n')

            # print(f':{data} :{status}')
    # print(leaf_types) 
    f_out.write(tup[2]) 
    f_out.close()

    # Do tag: Declare Branches
    tup = tag_value('src/events.h','Declare Branches')
    f_out = open('src/events.h','w')
    f_out.write(tup[0])
    f_out.write('    // List of branches\n')
    for data, status in zip(branch_names, keep):
        if status == 'keep':
            f_out.write(f'    TBranch        *{data};   //!\n')
        # don't write branches for TClonesArray or TObject instances
            # continue
        # elif status == 'keep':
        # else:
            # f_out.write(f'    TBranch        *b_{status};   //!\n')
        # print('this',data,status)
    f_out.write(tup[2])
    f_out.close()
    
    # Do tag: Set Branches
    tup = tag_value('src/events.cxx','Set Branches')
    f_out = open('src/events.cxx','w')
    f_out.write(tup[0])
    f_out.write("fChain->SetMakeClass(0); // note if is SetBranchAddress(1) then TObject's\n")
    f_out.write("                         // couldn't be read from the tree\n")
    # initialize the TClonesArrays
    # for key in clones:
        # print (best_matches_clones[key])
        # f_out.write(f'    {key:15s} = new TClonesArray("{best_matches_clones[key][0][0]}");\n');
    for key in classes:
        f_out.write(f'    {key:15s} = nullptr;\n');
    f_out.write('    // List of branches\n')
    for data, status in zip(set_branches, keep):
        if status == 'skip':
            continue
        elif status == 'keep':
            f_out.write(f'    fChain->SetBranchAddress("{data[0]}", {data[1]}, {data[2]});\n')
        else:
            if status in clones:
            # print(f'data({data}) status({status}) branch({branch})')
                f_out.write(f'    fChain->SetBranchAddress("{status}", &tca_{status});\n') #, &b_{status});\n')
            else:
                f_out.write(f'    fChain->SetBranchAddress("{status}", &{status});\n') #, &b_{status});\n')
            # f_out.write(f'    fChain->SetBranchAddress("{data[0]}", {data[1]}, {data[2]});\n')
    f_out.write(tup[2])
    f_out.close()

    # add the psueo-iterator code
    # set_tag_value('./src/events.cxx','Coda Functions',
    text = ''
    for status in keep: # this iterator is used to keep the same order as declarations
        if status in clones:
            _class = best_matches_clones[status][0][0]
            text += f'''
// Two  psuedo-iterator functions for TClonesArray({_class}) *tca_{status}
{_class}* events::get_{status}(int i) {{
    return ({_class}*) tca_{status}->UncheckedAt(i);
}};
'''
    set_tag_value('./src/events.cxx','Coda Functions',text) 
# swap_trees('alt')
import shutil
# shutil.copyfile('../pure/src/events.h','src/events.h')
# shutil.copyfile('../pure/src/events.cxx','src/events.cxx')

re_leaf_wTObjects('./','${HOME}/AN_common/include/TreeObj.h')
