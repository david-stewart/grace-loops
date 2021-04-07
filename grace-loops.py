#!/gpfs/loomis/apps/avx/software/Python/3.6.4-foss-2018a/bin/python3

'''
    This script setups directory to write loops to run over various local trees,
    here on GRACE, in parallel. 

    It may be run in three modes:
    A. -s, --setup
    B. -d, --hadd
    C. -a, --add-loop
    D. -c, --copy-loop

    These modes are as follow:
    ------------
    | A. setup |
    ------------
    Input:
       $ grace-loops.py -s PATH
       -or equivalent-
       $ grace-loops.py --setup PATH

       Uses options with defaults:
           --tree-name  (default: 'events') TTree name in ROOT files
           --n-files    (default: 10)       number of jobs per SLURM submission

       PATH may either be a directory containing the ROOT files, or a file containing that path

       Can change option defaults values ("new_tree_name" and "50" here):
       $ grace-loops.py -s PATH -t new_tree_name -n 50
        

    Output:
        Generates file tree:
            .
            ├-- Makefile
            ├-- sub-slurm (slurm submit scripts)
            |   |-- test_loop.sh
            ├-- out-slurm/ (dir of slurm output files)
            ├-- src
            |   ├-- events.h   (with data from ./class.{h,C})
            |   ├-- events.cxx (with data from ./class.{h,C})
            |   ├-- main.cxx 
            |   ├-- MemTimeProgression.cxx 
            |   ├-- MemTimeProgression.h 
            |   ├-- lists.h 
            |   ├-- lists.cxx
            |   |-- test_loop.cxx
            ├-- in-lists       
            |   ├-- list_0.list (each list with n_files)
            |   ├-- list_1.list
            |   |-- ...
            |   |-- list_all.list
            ├-- bin/
            ├-- obj/
            ├-- out-data
            |   |-- test_loop (output for loop "test_loop")
            |       |-- out-files
            ├-- _x.cxx -> ./src/test_loop.cxx (soft link)
            ├-- _m.cxx -> ./src/main.cxx      (soft link)
            |-- _h.h   -> ./src/events.h      (soft link)
            |-- fastjet-> ${FASTJET3}/tools/fastjet (soft link)

    ------------
    | B. hadd  |
    ------------
    To run use option: -d, --hadd

    Input:
       $ grace-loops.py -d LOOP_NAME
       -or equivalent-
       $ grace-loops.py --hadd  LOOP_NAME

    Output:
           Generates hadd file:
            ├-- out-data
            |   |-- LOOP_NAME
            |       ├-- hadd.root <- This file
            |       |-- out-files
            |           |-- ...
        
    ---------------
    | C. add-loop |
    ---------------
    To run use option: -a, --add-loop

    Input:
       $ grace-loops.py -a LOOP_NAME
       -or equivalent-
       $ grace-loops.py -a LOOP_NAME
       -or equivalent-
       $ grace-loops.py --add-loop  LOOP_NAME

    Output:
            Add/modifies files as follows:
            ├-- sub-slurm
            |   |-- LOOP_NAME.sh (new file)
            ├-- src
            |   ├-- events.h   (add friend function for function LOOP_NAME)
            |   ├-- main.cxx   (add call function for LOOP_NAME)
            |   |-- LOOP_NAME.cxx (skeleton code for new loop)
            ├-- Makefile (modified to compile with ./src/LOOP_NAME.cxx)
            ├-- _x.cxx -> ./src/LOOP_NAME.cxx (update soft link)
            |-- out-data
                |-- LOOP_NAME
                    |-- out-files


    ----------------
    | D. copy-loop |
    ----------------
    Input:
        $ grace-loops.py --copy-loop LOOP_FROM LOOP_TO

    Output:
        Will generate a new loop LOOP_TO just like $grace-loops -a LOOP_TO, but
        will copy over the code from ./src/LOOP_FROM.cxx and
        ./sub-slurm/LOOP_FROM.sh so that you can edit a new copy of the
        LOOP_FROM loop.
        

    --------------
    | User usage |
    --------------
        1. $ grace-loops.py --setup PATH
        2. edit _x.cxx -> src/test_loop.cxx
        3. $ sbatch sub-slurm/test_loop.sh
           generates files in out-data:
            ├-- out-data
            |   |-- test_loop
            |       |-- out-files
            |          ├-- test_loop-1.root
            |          ├-- test_loop-2.root
            |          |-- ...
        4. $ grace-loops.py hadd test_loop
        5. use ./out-data/test_loop/hadd.root
        6. For as many user loops as needed:
           *  $ grace-loops.py --add-loop LOOP_NAME
              -or-
              $ grace-loops.py --copy-loop LOOP_FROM LOOP_TO
           *  repeat steps 2.-5.
'''

from string import Template
import subprocess, os, argparse
from sys import argv
from glob import glob
from pathlib import Path
import re
import datetime

class file_times:
    def __init__(self, glob_string):
        # print(f'glob_string {glob_string}')
        self.youngest_time = -1
        self.files = []
        # self.times = []
        for F in glob(glob_string):
            # print(F)
            time = os.stat(F).st_mtime
            # print(time)
            # self.times.append(time)
            self.files.append((time,F))
            if time > self.youngest_time:
                self.youngest_time = time
        self.nfiles = len(self.files)
        self.files.sort(reverse=True)
        # print(self.nfiles, self.youngest_time)
    def fmt(time_f):
        time = datetime.datetime.fromtimestamp(time_f)
        return f'{time.year}-{"%02i"%time.month}-{"%02i"%time.day} {"%02i"%time.hour}:{"%02i"%time.minute}:{"%02i"%time.second}'

def gen_dict():
    ''' Make the input dictionary used in generating the input files '''
    return {

'Makefile' : 
r'''
# The libraries for pythia, fastjet, and root ( "LIB_TRI" )
# ccflg=-I${PYTHIA8}/include `${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags`
# ldflg=-L${PYTHIA8}/lib -lpythia8 `${FASTJET3}/fastjet-config --libs` `root-config --glibs`
 
# ccflg= `root-config --cflags`
ccflg=`${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags`
# ldflg= `root-config --glibs`
ldflg=`${FASTJET3}/fastjet-config --libs` `root-config --glibs`

# LIB_PYTH=-I${PYTHIA8}/include -L${PYTHIA8}/lib -lpythia8
LIB_FASTJET=`${FASTJET3}/fastjet-config --cxxflags --libs`
LIB_ROOT=`root-config --cflags --glibs`
LIB_DUO= ${LIB_ROOT} ${LIB_FASTJET}
# LIB_TRI=${LIB_PYTH} ${LIB_FASTJET} ${LIB_ROOT}

# compilation option
CC=g++
CFLAGS=-std=c++11 -O3 -Wno-deprecated
CFLAGS_CHECK=-std=c++11 -O0 -Wno-deprecated -g

bin/main: obj/events.o \
          obj/main.o   \
          obj/MemTimeProgression.o
	${CC} ${CFLAGS} -o $@ $^ ${LIB_DUO} 

obj/events.o: src/events.cxx src/events.h src/lists.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/main.o: src/main.cxx src/events.h src/MemTimeProgression.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@

obj/MemTimeProgression.o: src/MemTimeProgression.cxx src/MemTimeProgression.h
	${CC} ${CFLAGS} ${ccflg} -c $< -o $@
''',

'MemTimeProgression.cxx':
'''
#include "MemTimeProgression.h"

#include "stdlib.h"
#include "stdio.h"
#include <string>
#include <iostream>
#include "TString.h"

/* ClassImp(MemTimeProgression) */

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\\0';
    i = atoi(p);
    return i;
};

int getMemValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
};

MemTimeProgression::MemTimeProgression(int print_int=200) :
    nCalls   {0},
    mem0     {0},
    max_mem  {0},
    time0    {0.},
    watch    {},
    call_print_interval {print_int}
{
    watch.Start();
    call();
};

string MemTimeProgression::set_stats() {
    int    mem1  = getMemValue();
    if (mem1 > max_mem) max_mem = mem1;
    double time1 = watch.RealTime();
    watch.Continue();

    const char* pm_mem = (mem1>mem0) ? "+" : "-";
    stats = Form(" Finished %8lli calls | Time: %5.0f sec (+ %4.0f) | "
                 "Mem: %6.2f MB (%s%6.2f)",
                  nCalls, time1, time1-time0, mem1/1000., pm_mem, (mem1-mem0)/1000.);
    time0=time1;
    mem0=mem1;
    return stats;
};

bool MemTimeProgression::call(){
    ++nCalls;
    if (nCalls % call_print_interval == 0) {
        set_stats();
        cout << stats << endl;
        return true;
    } else {
        return false;
    }
};

string MemTimeProgression::set_get_stats() {
    cout << set_stats() << endl;
    return stats;
};
''',

'lists.h' :
'''#ifndef lists__h
#define lists__h

#include <vector>
using namespace std;
// need functions to read in the runId lists
vector<int> read_file_list(const char* list, int col=0, bool sort=false);

#endif
''',

'lists.cxx' :
'''#include "lists.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <iostream>

using namespace std;

vector<int> read_file_list(const char* file_name, int col, bool do_sort){
    ifstream f_in;
    f_in.open(file_name);
    vector<int> vec {};
    if (!f_in.is_open()) {
        cout << " fatal: couldn't open input file " << file_name;
        return vec;
    }
    string line;
    while (getline(f_in, line)){
        if (line.rfind("#",0) == 0) continue; // it begins with # and is a comment

        istringstream is{ line };
        int i_val;
        for (int i{0};i<=col;++i) {
            is >> i_val;
            /* cout << " i: " << i << " val " << i_val << endl; */

        }
        vec.push_back(i_val);
    }
    f_in.close();

    // sort the list
    if (do_sort) sort(vec.begin(), vec.end());
    return vec;
};
''',

'MemTimeProgression.h':
'''
#ifndef MemTimeProgression__h
#define MemTimeProgression__h

// determine how much memory has been used

#include "stdlib.h"
#include "stdio.h"
#include <string>
#include "TStopwatch.h"

int parseLine(char* line);
int getMemValue();
using namespace std;

struct MemTimeProgression {
    //--------------------------------
    // start at step=0
    // increment with each call()
    //  -> print to cout at call_print_interval
    //  return true is printing (to allow to print out to log)
    MemTimeProgression(int step_int);
    long long int nCalls;
    int    mem0;
    int    max_mem;
    double time0;
    TStopwatch watch;
    int  call_print_interval; 
    bool call(); // print every print_interval calls; return true if print
    string set_stats();
    string set_get_stats();
    string stats; // populate message
};

#endif
''',

'events.h' : Template( # template keys: 
                       # - $KEY_constexpr
                       # - $KEY_leaf_types
                       # - $KEY_branches
'''
#ifndef events_h
#define events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"
#include "TClonesArray.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include "MemTimeProgression.h"

class events {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   long long int nevents;
   ofstream &log;
   MemTimeProgression stats;

//:TAG START:  Array Sizes
$KEY_constexpr
//:TAG END: Array Sizes

//:TAG START: Leaf Types
$KEY_leaf_types
//:TAG END: Leaf Types

//:TAG START: Declare Branches
$KEY_branches
//:TAG END: Declare Branches

   events(ofstream& log, int n_events, TString inlist);
   virtual ~events();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // values for running the loop
   bool next();
   Long64_t nentries;
   Long64_t jentry;

   // Make friend functions to fill in the runs
   // i.e. friend void runLoop(events&, string);
   // TAG: start-friend-functions

};
#endif
'''),

'events.cxx' : Template( # template key:
                         # - KEY_SetBranchAddress
                         # - KEY_TreeName
'''
#define events_cxx
#include "events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


events::events(ofstream& _log, int _n_events, TString inp_file) :
    nevents{_n_events},
    log {_log}, 
    stats {1000},
    jentry {-1}
{
    TChain* tree  = new TChain("$KEY_TreeName");
    if (inp_file.EndsWith(".root")) {
            cout << " Adding input file: " << inp_file << endl;
        tree->Add(inp_file.Data());
    } else if (inp_file.EndsWith(".list")) {
        string line;
        ifstream list;
        list.open(inp_file.Data());
        while (getline(list, line)){
            cout << " Adding input file: " << line << endl;
            tree->Add(line.c_str());
        }
        list.close();
    }
    /* cout << " has " << tree->GetEntries() << " in tree" << endl; */
    /* cout << " has " << tree->GetEntriesFast() << " in tree" << endl; */
    Init(tree);
    nentries = tree->GetEntries();
    if (nevents == -1) nevents = nentries;
}
events::~events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t events::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t events::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //:TAG START: Set Branches
   fChain->SetMakeClass(1);

$KEY_SetBranchAddress
   //:TAG END: Set Branches
   Notify();
}

Bool_t events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void events::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t events::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

bool events::next() {
    jentry++;
    if (jentry >= nevents) {
        stats.set_get_stats();
        log  << " Final stats: " << stats.stats << endl; 
        return false;
    }
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
        cout << " Breaking out of loop at jentry on failure to read ientry: " << jentry << endl;
        stats.set_get_stats();
        log  << " Final stats: " << stats.stats << endl; 
        return false;
    }
    fChain->GetEntry(jentry);
    if (stats.call()) log << stats.stats << endl;
    return true;
}

//:TAG START: Coda Functions
//:TAG END: Coda Functions

'''),

'loop.sh': Template( # keys:
                     # - KEY_loop_name
                     # - KEY_n_lists
'''#!/usr/bin/bash
#SBATCH  --job-name=$KEY_loop_name
#SBATCH  --array=0-$KEY_n_lists
#SBATCH  -o out-slurm/%a

#SBATCH   --partition=day
#SBATCH   --time=00:10:00
#SBATCH   --mem=1G

n_events=-1
which_loop="$KEY_loop_name"
dir="./out-data/$KEY_loop_name/out-files"
inp_list="./in-lists/list_$${SLURM_ARRAY_TASK_ID}.list"
name="$${which_loop}"

./bin/main $${n_events} $${inp_list} $${name} $${dir}/$${name}_$${SLURM_ARRAY_TASK_ID}
'''
),

'main.cxx' :
'''
#include "src/events.h"
#include <sstream>
#include <fstream>
#include "TFile.h"
// #include <algorithm>
using namespace std;

int main(int nargs, char** argv) {
    /*
     * arguments:
     *   1: number of events
     *   2: input root file list name
     *   3: which program to run
     *   4: output base-name
     *   5: optional input
     */

    int n_events          {  (nargs>1) ? atoi(argv[1]) : 100 };
    string inp_list       {  (nargs>2) ? argv[2] : "in-lists/list_all.list" };
    string which_loop     {  (nargs>3) ? argv[3] : "test_loop" };
    string o_name_tag     {  (nargs>4) ? argv[4] : "test_loop" };

    ostringstream collect;
    for (int i{5};i<nargs;++i) {
        string arg {argv[i]};
        collect << arg << " ";
    }

    ofstream log;
    log.open((o_name_tag + ".log").c_str());
    log << "Starting output."  << endl
        << "Command line input:" << endl;
    for (int i{0}; i<nargs; ++i) log << "arg("<<i<<")  " << argv[i] << endl;
    log << endl << endl;

    events my_events{log, n_events, inp_list};

    
    TFile fout  { (o_name_tag+".root").c_str(), "recreate" };

    // run the loop
    cout << " Looking for function: " << which_loop << endl;
    if (which_loop == "empty-loop") {
        // TAG: empty-loop
    } else {
        cout << " -  Fatal error: couldn't find loop \\"" << which_loop << "\\"" << endl
            << "     ->  Terminating program." << endl;
    }
    fout.Close();
    log.close();
};
''',

'loop.cxx' : Template( # key: KEY_loop_name
'''
#include "events.h"
#include <sstream>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

using namespace std;
void $KEY_loop_name(events& dat, string _options) {
    // Start of 
    cout << " Running fn \\"$KEY_loop_name\\"" << endl;
    istringstream options ( _options );
    int n_options = 0;
    string arg;
    // options >> arg;
    while (options >> arg) {
        cout    << " Option " << n_options << ":  " << arg << endl;
        dat.log << " Option " << n_options << ":  " << arg << endl;
        ++n_options;
    }

    // Histogram declarations here:
    // TH1D hg {"hg", "a;b;c", 10, 0., 1.};


    // Run loop here:
    while (dat.next()) {
        // cout << " Doing nothing in event " << dat.jentry << endl;
        // dat.log  << " Doing nothing in event " << dat.jentry << endl;
        // hg.Fill( x );
    }

    // Write histograms here
    // hg.Write();


    // Wrap-up work here:
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
'''),

}

def fix_double_branches(file_name):
    '''Under the certain conditions, ROOT's TTree->MakeClass("name") will
       make multiple TBranchs of the same name. This is cause c++ to not 
       compile the code as desired.

       Input:
           file name
       Process:
           Check all the TBranch lines for double names and rename duplicates to unique values
       Output:
           Number of duplicates found and fixed       
   '''
    if not os.path.isfile(file_name):
        print(f'fatal: cannot find input file {file_name} in function "fix_double_branches"')
        exit()

    collect_lines = []
    n_subs = 0;
    name_set = set();
    re_name = None
    if (file_name[-2:] == '.h'):
        re_name = re.compile(r'(\s*TBranch\W+\*)(\w+)(;.*)')
    elif (file_name[-4:] == '.cxx' or file_name[-2:] == '.C'):
        re_name = re.compile(r'(\s*fChain->SetBranchAddress.*\&)(\w+)(.*)')
    else:
        print('fatal: argument to fix_double_branches doesn\'t end in .h .cxx or .C')
        exit()

    for line in open(file_name,'r').readlines():
        line = line.rstrip();
        match = re_name.match(line)
        if match:
            name = match.group(2)
            if name in name_set:
                n_subs += 1
                i = 0
                prop_name = f"{name}_{i}"
                while prop_name in name_set:
                    i += 1
                    prop_name = f"{name}_{i}"
                name = prop_name
            collect_lines.append(f'{match.group(1)}{name}{match.group(3)}')
            name_set.add(name)    
        else:
            collect_lines.append(line)
    if n_subs:
        with open(file_name,'w') as f_out:
            f_out.write('\n'.join(collect_lines))
    return n_subs

# if __name__ == '__main__':
    # n_subs = fix_double_branches("src/events.cxx")     
    # print(f'n_subs: {n_subs}')

#   ------------------------------
#   | program module: "add-loop" |
#   ------------------------------
#    note: this is called by "setup" to make the initial loop "test_loop"
#          and therefore must be defined first
def add_loop(loop_name, templates=None):
    for F in ('Makefile','src/main.cxx','src/events.h'):
        if not os.path.isfile(F):
            print(f'fatal: required file "{F}" not present for add-loop')
            exit(2)

    if os.path.isfile(f'src/{loop_name}.cxx'):
        print(f'fatal: loop file  src/{loop_name}.cxx is already present.')
        exit(2)

    if not templates:
        templates = gen_dict()

    keys = {'KEY_loop_name' : loop_name }
    imax_list = -1;
    while os.path.isfile(f'in-lists/list_{imax_list+1}.list'):
        imax_list += 1
    keys['KEY_n_lists'] = imax_list
    
    with open(f'sub-slurm/{loop_name}.sh','w') as fout:
        fout.write(templates['loop.sh'].substitute(**keys))
    os.symlink(f'sub-slurm/{loop_name}.sh',f'{loop_name}.sh')

    

    # modify src/events.h in place
    text = Path('src/events.h').read_text()
    parts = re.search('(.*TAG: start-friend-functions\n)(.*)',text,re.DOTALL)
    with open('src/events.h','w') as fout:
        fout.write(
        parts.group(1) + f"   friend void {loop_name}(events&, string);\n" + parts.group(2) )

    # modify src/main.cxx in place
    text = Path('src/main.cxx').read_text()
    parts = re.search(r'(.*\n\s*string which_loop[^"]*")\w*'
        +r'(".*\n\s*string o_name_tag[^"]*")\w*'
        +r'(.*// TAG: empty-loop\w*\n)(.*)', text,re.DOTALL)
    with open('src/main.cxx','w') as fout:
        fout.write(parts.group(1))
        fout.write(loop_name)
        fout.write(parts.group(2))
        fout.write(loop_name)
        fout.write(parts.group(3))
        fout.write(f'''   }} else if (which_loop == "{loop_name}") {{
        {loop_name}(my_events, collect.str());\n''')
        fout.write(parts.group(4))

    # write the new src/loop_name.cxx file
    with open(f'src/{loop_name}.cxx','w') as fout:
        fout.write(templates['loop.cxx'].substitute(**keys))

    # modify the Makefile in place:
    text = Path('Makefile').read_text()
    parts = re.search(r'(.*\nbin/main:.*?\n.*?)(\n\t\${CC}.*)',text,re.DOTALL)
    with open('Makefile','w') as fout:
        fout.write('\n\n')        
        fout.write(parts.group(1))
        fout.write(f' \\\n          obj/{loop_name}.o')
        fout.write(parts.group(2))
        fout.write(f'\nobj/{loop_name}.o: src/{loop_name}.cxx src/events.h\n')
        fout.write( '\t${CC} ${CFLAGS} ${ccflg} -c $< -o $@\n')

    # update the _x.cxx link
    subprocess.Popen(['ln','-fs', f'src/{loop_name}.cxx', '_x.cxx'])

    # add output directories
    for D in (f'out-data/{loop_name}',f'out-data/{loop_name}/out-files'):
        if os.path.isdir(D):
            print(f'Warning: directory "{D}" already exists')
        else:
            os.mkdir(D)

def copy_loop(loop_from, loop_to, templates=None):
    # Check that the loop_from exists
    if not os.path.isfile(f'src/{loop_from}.cxx'):
        print(f'fatal: loop from file src/{loop_from}.cxx does not exist')
        exit()

    print(f'Warning: all instances of word "{loop_from}" will \n'
           'be replaced in "{loop_to}" in files ./src/{loop_to} \n'
           'and ./sub-slurm/{loop_to}, regardless of correctness.\n'
           '-> Make sure "{loop_from}" is unique in these files.')

    #FIXME
    add_loop(loop_to)

    # Now fix the existing files:
    # sub-slurm/{loop_from.sh}
    for file in ('src/NAME.cxx', 'sub-slurm/NAME.sh'):
        txt = Path(file.replace('NAME',loop_from)).read_text()
        txt.replace(loop_from,loop_to)
        with open(file.replace('NAME',loop_to),'w') as f_out:
            f_out.write(txt.replace(loop_from,loop_to))

#   ---------------------------
#   | program module: "setup" |
#   ---------------------------
def setup(in_path, tree_name='events', max_nfiles=10, array_files=-1):
    '''See notes at start of file'''
    #---------------------------------
    #| Generate file tree structure  |
    #---------------------------------
    if os.path.isdir("src"):
        print("Directory \"src\" is already present. Terminating grace_loops.py setup")
        exit()

    with open("in-path",'w') as fout:
        fout.write(f'# setup with input path:\n{in_path}')

    for D in ('sub-slurm',
              'out-slurm',
              'src',
              'in-lists',
              'bin',
              'obj',
              'out-data'):
        if not os.path.isdir(D):
            os.mkdir(D)

    #------------------------------
    #| Write the input file lists |
    #------------------------------

    # set max_nfiles to be either:
    #   the number of input files
    #   the input files in a list
    if not os.path.isdir(in_path) and os.path.isfile(in_path):
        for lines in open(in_path,'r').readlines():
            L = lines.strip()
            if len(L) == 0 or L[0] == '#':
                continue
            print(f'Found input directory {L} in file {in_path}')
            in_path = L
            break
        # see if in_path is a file that contains the path
    if not os.path.isdir(in_path):
        exit(f'fatal: input directory {in_path} not found')

    in_files = glob(f'{in_path}/*.root')
    with open ('in-lists/list_all.list','w') as f_out:
        for f in in_files:
            f_out.write(f'{f}\n')

    n_lists    = 0
    n_per_file = 0
    n_plus_one = 0

    n_files = len(in_files)
    if max_nfiles >= n_files:
        n_lists = n_files
        n_per_file = 1
    else:
        n_lists = max_nfiles
        n_per_file = n_files//max_nfiles
        n_plus_one = n_files%max_nfiles

    index = 0
    for i_list in range(n_lists):
        with open(f'in-lists/list_{i_list}.list','w') as fout:
            n = n_per_file
            if i_list < n_plus_one:
                n += 1
            for i in range(index, index+n):
                fout.write(f'{in_files[i]}\n')
            index += n
        

    # if len(in_files) == 0:
    #     exit(f'fatal: no input *.root files in {in_path}')

    # file_cnt = 0
    # entry_cnt = 0
    # fout_all = open(f'in-lists/list_all.list','w')
    # for name_file in in_files:
    #     if (entry_cnt == 0):
    #         fout = open(f'in-lists/list_{file_cnt}.list','w')
    #         file_cnt += 1
    #     fout.write(f'{name_file}\n')
    #     fout_all.write(f'{name_file}\n')
    #     entry_cnt += 1
    #     if entry_cnt == max_nfiles:
    #         entry_cnt = 0
    # fout_all.close()


    #------------------------------------------------
    #| Make events.{C.h} |
    #------------------------------------------------
    if array_files == -1 :
        print(' Using TChain->MakeClass(events) to make file structure reading \n',
             f' all {len(in_files)} input files. Array sizes will match maximum size\n',
              ' of any file.\n\n',
              ' If this takes too long (from many input files), rerun with --n-array-files #\n',
              ' The program will find the array maximums from only this number files, and then\n',
              ' pick a conservatively larger value. Still, be careful of seg faults if these\n',
              ' are not large enough.')
    else:        
        print(' Using TChain->MakeClass(events) to make file structure reading \n',
             f' first {array_files} input files. The program will find the array \n',
              ' the maximums from only these files, and then pick a conservatively \n',
              ' larger value. Still, be careful of seg faults if these\n',
              ' are not large enough.')

    # make the class.h and class.C from root
    if not os.path.isfile('events.h'):
        _out = subprocess.Popen([f'''
root -l -b <<EOF
    TChain *chain = new TChain("{tree_name}");
    ifstream fin;
    fin.open("in-lists/list_all.list");
    string file_path;
    int n_files = 0;
    int to_read = {array_files};
    while (fin >> file_path) {{
        chain->Add(file_path.c_str());
        if ( to_read != -1 && n_files >= to_read) break;
        ++n_files;
    }}
    chain->MakeClass("events")
EOF
    '''],shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)

        stdout,stderr = _out.communicate()
        try:
            stdout = stdout.decode('utf-8')
            stderr = stderr.decode('utf-8')
        except:
            pass
        print('\n\n--begin: ROOT messages from MakeClass("events")--')
        print(stdout)
        print('--end: ROOT messages from MakeClass("events")--\n\n')

        for suffix in ('C','h'):
            if not os.path.isfile(f'events.{suffix}'):
                # print(f'fatal: file events.{suffix} not generated')
                print(f'out msg:\n {stdout}')
                print(f'err msg:\n {stderr}')

        if not os.path.isfile('events.h'):
            exit('fatal: events.h not generated in time for script.\n'
                 ' please run script again.')


    
    #------------------------------------------------
    #| Get the required code snippets from events.h |
    #------------------------------------------------
    text = Path('events.h').read_text()

    # find the constexpr lines like:
    #     static constexpr Int_t kMaxtrack = 83;
    code_snippets = dict()
    static_constexpr = re.compile('(   static constexpr Int_t \w+ = )(\d*);')
    KEY_constexpr = []
    first = True
    mult_up = (array_files != -1)
    for x in static_constexpr.findall(text):
        if first and mult_up:
            first = False
            print('Finding max values for arrays and setting them to 4x(max(max_val, 20))')
        val = int(x[1])
        if mult_up:
            val0 = val
            if val < 20:
                val = 20
            val *= 4
            print(f'{x[0]}   {val0} -> {val}')
        KEY_constexpr.append(f"{x[0]}{val};")
    code_snippets['KEY_constexpr'] = "// Fixed size dimensions of " +\
            "array or collections " +\
            "stored in the TTree if any.\n"+"\n".join(KEY_constexpr)

    # // Get the leaf types
    leaf_types = re.search(
        '(  // Declaration of leaf types.*)(\n   // List of branches)',
        text, re.DOTALL)
    code_snippets['KEY_leaf_types'] = leaf_types.group(1)
    # print(KEY_leaf_types)

    # Get list of branches
    KEY_branches = "   // List of branches\n"
    KEY_branches += "\n".join([ x.group() for x in re.finditer(
        '   TBranch .*',text)])
    code_snippets['KEY_branches'] = KEY_branches
    # print(KEY_branches)

    # Get SetBranchAddress code
    KEY_SetBranchAddress = "\n".join([x.group() for x in re.finditer(
        '   fChain->SetBranchAddress.*',text)])
    code_snippets['KEY_SetBranchAddress'] = KEY_SetBranchAddress
    code_snippets['KEY_TreeName'] = tree_name
    # print(KEY_SetBranchAddress)

    #-------------------------
    #| Write the input files |
    #-------------------------
    templates = gen_dict()

    with open('Makefile','w') as fout:
        fout.write(templates['Makefile'])
    
    with open('./src/events.h','w') as fout:
        fout.write(templates['events.h'].substitute(**code_snippets))

    with open('./src/events.cxx','w') as fout:
        fout.write(templates['events.cxx'].substitute(**code_snippets))

    with open('./src/main.cxx','w') as fout:
        fout.write(templates['main.cxx'])

    for suffix in ('cxx','h'):
        with open(f'./src/MemTimeProgression.{suffix}','w') as fout:
            fout.write(templates[f'MemTimeProgression.{suffix}'])
        with open(f'./src/lists.{suffix}','w') as fout:
            fout.write(templates[f'lists.{suffix}'])

    #-------------------
    #| make soft links |
    #-------------------
    for pair in (('src/events.h','_h.h'),
                ('src/main.cxx','_m.cxx'),
                ('${FASTJET3}/tools/fastjet','fastjet') ):
        subprocess.Popen(['ln','-fs', pair[0], pair[1]])

    #-------------------------
    #| add the test_loop     |
    #-------------------------
    add_loop ("test_loop", templates)

    #----------------------------------
    #| Fix potential double  branches |
    #----------------------------------
    nsubs = fix_double_branches('src/events.h')
    if nsubs > 0:
        fix_double_branches('src/events.cxx')

def hadd(name=None, o_tag=None):
    if not name:
        text = Path('Makefile').read_text()
        last_dot_o = text.rfind('.o')
        last_obj   = text.rfind('obj/')
        name = text[last_obj+4:last_dot_o]
        print(f'Attemping to run hadd in output "out-data/{name}/"')
    else:
        print(f'Running hadd in output "out-data/{name}/"')

    if o_tag:
        o_tag = '_'+o_tag
    else:
       o_tag = '' 

    hadd_files = glob(f'out-data/{name}/out-files/*.root')        
    if not (hadd_files):
        print (f'fatal: no out-data/{name}/out-files/*.root files to hadd')
        exit(1)


    if os.path.isfile(f'out-data/{name}/hadd{o_tag}.root'):
        print(f'warning: writing over out-data/{name}/hadd{o_tag}.root')

    _out = subprocess.Popen(f'hadd -fkT out-data/{name}/hadd{o_tag}.root {" ".join(hadd_files)}',shell=True)
    stdout,stderr = _out.communicate()
    try:
        stdout = stdout.decode('utf-8')
        print(f'out msg:\n {stdout}')

        stderr = stderr.decode('utf-8')
        print(f'err msg:\n {stderr}')
    except:
        pass

    tree_obj = []
    try:
        _out = subprocess.Popen(f'''root -l {hadd_files[0]}<<EOF  > __temp_types
.ls
EOF
''', shell=True)
        for line in open('__temp_types','r').readlines():
            if 'KEY: TTree' in line:
                tree_obj.append(line.split(';')[0].split()[-1])
    except:
        pass 
#-------------------------------
    for tree in tree_obj:
        _out = subprocess.Popen(
            f'rootcp {hadd_files[0]}:{tree} out-data/{name}/hadd{o_tag}.root',shell=True);
    # except:
        # pass

    print('suggested:')
    print(f'scp grace:{os.getcwd()}/out-data/{name}/hadd{o_tag}.root .') 


#-------------------------
#| start main program    |
#-------------------------
def parse_args():
    parser = argparse.ArgumentParser(
            description="Required: run in one of 4 selected modes: [-s | -a | -d | -c]")
    parser_mode = parser.add_mutually_exclusive_group()
    # one of three required modes
    parser_mode.add_argument('-s','--setup',  help='input path or file-containing-path; '
            ' uses TREE_NAME & N_FILES')
    parser_mode.add_argument('-a','--add-loop',  help='add a loop')
    parser_mode.add_argument('-d','--hadd',      help='hadd an existing loop')
    parser_mode.add_argument('-c','--copy-loop', 
    help='copy an existing loop; requires name from and name to', nargs=2)

    # individual arguments
    parser.add_argument('-t','--tree-name', 
            help='Name of TTree in ROOT files. Default:"events"', type=str, default="events")
    parser.add_argument('-n','--n-files',   
            help="Divide slurm submissions into this many jobs. Default:10", type=int, default=10)
    parser.add_argument('--tag',            help='tag name for hadd')
    parser.add_argument('--n-array-files',   
            help="Number of files read to TChain for finding max array sizes. Default:20. Use -1 for all", type=int, default=15)
    parser.add_argument('-p','--print-path', help='output for copy-paste; if want all use "."')
    args = parser.parse_args()

    if not (args.setup or args.add_loop or args.hadd or args.copy_loop or args.print_path):
        parser.print_help()
        exit()

    # print(args)
    return args

def print_path(which=""):
    if which[0] == '+':
        by_time = False
        which = which[1:]
    else:
        by_time = True
        

    pwd = os.getcwd()

    if not which == '.':
        ofile = f'{pwd}/out-data/{which}/hadd.root'
        if not os.path.isfile(f'ofile'):
            print(f'warning: file {ofile} is not present')
            print_path('.')
        else:
            print('ofile')
        return

    if not by_time:
        for ofile in glob(f"{pwd}/out-data/*/hadd*.root"):
            print('  ',ofile)
        return

    # make a sorted time list of all out-data files
    data = dict()
    for D in glob(f"{pwd}/out-data/*"):
        if os.path.isdir(D):
            _hadd_data = file_times(f'{D}/*.root')
            _root_data = file_times(f'{D}/out-files/*.root')
            if _hadd_data.nfiles and _root_data.nfiles:
                time = min(_hadd_data.youngest_time, _root_data.youngest_time)
                data[(time,D)] = (_hadd_data, _root_data)
            elif _hadd_data.nfiles:
                time = _hadd_data.youngest_time
                data[(time,D)] = (_hadd_data, None)
            elif _root_data.nfiles:
                time = _root_data.youngest_time
                data[(time,D)] = (None, _root_data)
    for key in sorted(data,reverse=True):
        E = data[key]
        if E[0] and E[1]:
            print( f'{file_times.fmt(E[1].youngest_time)}  |{E[1].files[0][1].split("/")[-3]}|  : {E[1].nfiles} input *.root files')
            for F in E[0].files:
                time = F[0]
                name = F[1]
                lead = "  !->" if time<E[1].youngest_time else "   ->"
                print( f'{lead} {file_times.fmt(time)} {name}')
        elif E[0]:
            print( f'{file_times.fmt(E[0].youngest_time)}  |{E[0].files[0][1].split("/")[-2]}|  : 0 *.root files')
            # print( f'{file_times.fmt(E[0].youngest_time)}  {os.path.dirname(E[0].files[0][1])}/<no-files>.root')
            for F in E[0].files:
                print( f'   -> {file_times.fmt(F[0])} {F[1]}')
        elif E[1]:
            print( f'{file_times.fmt(E[1].youngest_time)}  |{E[1].files[0][1].split("/")[-3]}|  : {E[1].nfiles} *.root files')
        # print(key[1])
    # print out the output

    # else:
    #     index = dict()
    #     index_y_file = dict()
    #     for ofile in o_files:
    #         time = os.stat(ofile).st_mtime
    #         # get the oldest sub-file, too
    #         youngest_file = -1.
    #         for F in glob(f'{os.path.dirname(ofile)}/*.root'):
    #             if youngest_file < os.stat(F).st_mtime:
    #                 youngest_file = os.stat(F).st_mtime
    #         index_y_file[youngest_file] = 
    #         index[os.stat(ofile).st_mtime] = ofile
    #     K = sorted(index.keys(),reverse=True)
    #     for k in K:
    #         time = datetime.datetime.fromtimestamp(k)
    #         print (f'{time.year}-{time.month}-{time.day} {time.hour}:{time.minute}:{"%02i"%time.second} : {index[k]}')
    #         # print (datetime.datetime.fromtimestamp(k))
    #         # print(k)
if __name__ == "__main__":
    args = parse_args()
    # print(args)
    if args.setup:
        setup(args.setup, args.tree_name, args.n_files,args.n_array_files)
    elif args.hadd:
        # print('system args: ', args.hadd, '  next: ', args.tag)
        # exit('done')
        hadd(args.hadd, args.tag)
    elif args.add_loop:
        add_loop(args.add_loop) 
    elif args.copy_loop:
        copy_loop(args.copy_loop[0], args.copy_loop[1]) 
    elif args.print_path:
        print_path(args.print_path)
    exit(0)
