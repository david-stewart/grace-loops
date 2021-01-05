#!/gpfs/loomis/apps/avx/software/Python/3.6.4-foss-2018a/bin/python3

'''
    This script setups directory to write loops to run over various local trees,
    here on GRACE, in parallel. 

    It may be run in three modes:
    A. setup
    B. hadd
    C. add loop

    These modes are as follow:
    ------------
    | A. setup |
    ------------
    Input:
       $ grace-loops.py setup <inp-root-files-path> <tree-name=events> <n_files per input=10>
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
    Input:
       $ grace-loops.py hadd <loop-name=last-loop-name-generate>
    Output:
           Generates hadd file:
            ├-- out-data
            |   |-- test_loop
            |       ├-- hadd.root <- This file
            |       |-- out-files
            |           |-- ...
        
        Run: `MakeClass()` on one of input root files trees:
    
    ---------------
    | C. add-loop |
    ---------------
    Input:
       $ grace-loops.py add-loop <loop-name>
    Output:
            Add/modifies files as follows:
            ├-- sub-slurm
            |   |-- <loop-name>.sh (new file)
            ├-- src
            |   ├-- events.h   (add friend function for function <loop-name>)
            |   ├-- main.cxx   (add call function for <loop-name>)
            |   |-- <loop-name>.cxx (skeleton code for new loop)
            ├-- Makefile (modified to compile with ./src/<loop-name>.cxx)
            ├-- _x.cxx -> ./src/<loop-name>.cxx (update soft link)
            |-- out-data
                |-- <loop-name>
                    |-- out-files


    --------------
    | User usage |
    --------------
        1. $ grace-loops.py <path-to-input> <tree-name=events> <files-per-input=10>
        2. edit _x.cxx
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
           1. $ grace-loops.py add-loop <loop-name>
           2. repeat steps 2.-5.

'''

from string import Template
import subprocess, os
from sys import argv
from glob import glob
from pathlib import Path
import re

def gen_dict():
    ''' Make the input dictionary used in generated the input files '''
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
                  nCalls, time1, time1-time0, mem0/1000., pm_mem, (mem1-mem0)/1000.);
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

$KEY_constexpr

$KEY_leaf_types

$KEY_branches

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
    stats {1000} 
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
    jentry = 0;
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
   fChain->SetMakeClass(1);

$KEY_SetBranchAddress
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
    jentry++;
    return true;
}
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
    string inp_list       {  (nargs>2) ? argv[2] : "in-lists/list_0.list" };
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
        cout    << " Option " << n_options++ << ":  " << arg << endl;
        dat.log << " Option " << n_options++ << ":  " << arg << endl;
    }

    // Histogram declarations here:


    // Run loop here:
    while (dat.next()) {
        cout << " Doing nothing in event " << dat.jentry << endl;
        dat.log  << " Doing nothing in event " << dat.jentry << endl;
    }


    // Wrap-up work here:
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;
}
'''),

}

#   ------------------------------
#   | program module: "add-loop" |
#   ------------------------------
#    note: this is called by "setup" to make the initial loop "test_loop"
#          and therefore must be defined first
def add_loop(loop_name, templates=None):
    if not templates:
        templates = gen_dict()

    keys = {'KEY_loop_name' : loop_name }
    keys['KEY_n_lists'] = len(glob('in-lists/*.list'))-1
    # templates['KEY_loop_name'] = loop_name 
    
    with open(f'sub-slurm/{loop_name}.sh','w') as fout:
        fout.write(templates['loop.sh'].substitute(**keys))

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

#   ---------------------------
#   | program module: "setup" |
#   ---------------------------
def setup(in_path, tree_name='events', max_nfiles=10):
    '''See notes at start of file'''

    #------------------------------------------------
    #| Find input *.root file and make events.{C.h} |
    #------------------------------------------------
    if not os.path.isdir(in_path):
        exit(f'fatal: input directory {in_path} not found')

    in_files = glob(f'{in_path}/*.root')
    if len(in_files) == 0:
        exit(f'fatal: no input *.root files in {in_path}')
    # for f in in_files:
    #     print(f)

    # make the class.h and class.C from root
    if not os.path.isfile('events.h'):
        _out = subprocess.Popen([f'''root -l -b {in_files[0]}<<EOF
    {tree_name}->MakeClass("events")
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

        for suffix in ('C','h'):
            if not os.path.isfile(f'events.{suffix}'):
                # print(f'fatal: file events.{suffix} not generated')
                print(f'out msg:\n {stdout}')
                print(f'err msg:\n {stderr}')

        if not os.path.isfile('events.h'):
            exit('fatal: events.h not generated in time for script.\n'
                 ' please run script again.')

    #----------------------------
    #| Generate tree structure  |
    #----------------------------
    for D in ('sub-slurm',
              'out-slurm',
              'src',
              'in-lists',
              'bin','obj',
              'out-data'):
        if not os.path.isdir(D):
            os.mkdir(D)

    #------------------------------
    #| Write the input file lists |
    #------------------------------
    file_cnt = 0
    entry_cnt = 0
    for name_file in in_files:
        if (entry_cnt == 0):
            fout = open(f'in-lists/list_{file_cnt}.list','w')
            file_cnt += 1
        fout.write(f'{name_file}\n')
        entry_cnt += 1
        if entry_cnt == max_nfiles:
            entry_cnt = 0

    #------------------------------------------------
    #| Get the required code snippets from events.h |
    #------------------------------------------------
    text = Path('events.h').read_text()

    # find the constexpr lines like:
    #     static constexpr Int_t kMaxtrack = 83;
    code_snippets = dict()
    static_constexpr = re.compile('(   static constexpr Int_t \w+ = )(\d*);')
    KEY_constexpr = []
    for x in static_constexpr.findall(text):
        print(x)
        val = int(x[1])
        if val < 20:
            val = 20
        val *= 4
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

    # try:    
    with open('./src/main.cxx','w') as fout:
        fout.write(templates['main.cxx'])
    # except:
        # os.remove('./src/main.cxx')
        # with open('./src/main.cxx','w') as fout:
            # fout.write(templates['main.cxx'])

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

def hadd(name=None):
    if not name:
        text = Path('Makefile').read_text()
        last_dot_o = text.rfind('.o')
        last_obj   = text.rfind('obj/')
        name = text[last_obj+4:last_dot_o]
        print(f'Attemping to run hadd in output "out-data/{name}/"')
    else:
        print(f'Running hadd in output "out-data/{name}/"')

    hadd_files = glob(f'out-data/{name}/out-files/*.root')        
    if not (hadd_files):
        print (f'fatal: no out-data/{name}/out-files/*.root files to hadd')
        exit(1)

    if os.path.isfile(f'out-data/{name}/hadd.root'):
        print(f'warning: writing over out-data/{name}/hadd.root')

    _out = subprocess.Popen(f'hadd -fk out-data/{name}/hadd.root {" ".join(hadd_files)}',shell=True)
    stdout,stderr = _out.communicate()
    try:
        stdout = stdout.decode('utf-8')
        print(f'out msg:\n {stdout}')

        stderr = stderr.decode('utf-8')
        print(f'err msg:\n {stderr}')
    except:
        pass


#-------------------------
#| start main program    |
#-------------------------
    
#----tempary for testing
print(f'argv:{argv}')
if len(argv) == 1:
    if not os.path.isfile('bin/main'):
        argv.append('setup')
        argv.append('/gpfs/loomis/scratch60/caines/djs232/AN/pAu_JetEmb/TrackEfficiency')
        argv.append('primary')
        argv.append(20)
    else:
        hadd()
        exit(0)

if argv[1] == 'add-loop':
    for F in ('Makefile','src/main.cxx','src/events.h'):
        if not os.path.isfile(F):
            print(f'fatal: required file "{F}" not present for add-loop')
            exit(2)
    if len(argv)<3:
        print('fatal: to run add-loop a third argument is required')
    loop_name = argv[2]
    if os.path.isfile(f'src/{loop_name}'):
        print(f'fatal: loop file  src/{loop_name}.cxx is already present.')
        exit(2)
    add_loop(loop_name)
elif argv[1] == 'hadd':    
    if len(argv) > 2:
        hadd(argv[2])
    else:
        hadd()
elif argv[1] == 'setup':
    if len(argv) == 3:
        setup(argv[2])
    elif len(argv) == 4:
        setup(argv[2],argv[3])
    elif len(argv) > 4:
        setup(argv[2],argv[3],int(argv[4]))
else:
    print(''' Script usage: There are options:
$ grace-loops.py setup <inp-root-files-path> <tree-name=events> <n_files per input=10>
$ grace-loops.py hadd <loop-name>
$ grace-loops.py add-loop <loop-name>

Also see documentation in start of script.
''')
    exit(0)
