//Dane R Taylor May 13, 2015
//
// This code computes the full WTM map for a given network and threshold T
//
// input: -i linklist.tsv = contains the edge list with node indices 0,1,2,3,...
//        -t T = this is the threshold value
//
// output -o activation_times.tsv = the matrix of activation times, where row j corresponds to a contagion initialized at node j
//
//
//    ./compute_WTM_map_v3 -t 0.1 -i LinkList.tsv -o activation_times_v3.tsv

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

char *infile   = NULL;
char *outfile  = NULL;
char *Tvalue = NULL;




double strict_str2double(char* str)
{
    char* endptr;
    double value = strtod(str, &endptr);
    if (*endptr) return 0;
    return value;
}



// output how to use this file
void usage(char *prog_name, const char *more) {
    cerr << more << endl;
    cerr << "usage: " << prog_name << " -t T -i input_file -j input_node_file -o outfile" << endl << endl;
    cerr << "read the graph and convert it to binary format." << endl;
    exit(0);
}

// save inputs to variables
void parse_args(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'i':
                    if (i==argc-1)
                        usage(argv[0], "Infile missing\n");
                    infile = argv[i+1];
                    i++;
                    break;
                case 't':
                    if (i==argc-1)
                        usage(argv[0], "T value missing\n");
                    Tvalue = argv[i+1];
                    i++;
                    break;
                case 'o':
                    if (i==argc-1)
                        usage(argv[0], "Outfile missing\n");
                    outfile = argv[i+1];
                    i++;
                    break;
                default:
                    usage(argv[0], "Unknown option\n");
            }
        } else {
            usage(argv[0], "More than one filename\n");
        }
    }
    if (infile==NULL || outfile==NULL)
        usage(argv[0], "In or outfile missing\n");
}

// main to run
int main(int argc, char **argv) {
    clock_t t;
    
    parse_args(argc, argv);
    
    double T = strict_str2double(Tvalue);
  
    
    ////////////////////////////////////////////////////////////////////////////////
    // load link list and generate network
    ////////////////////////////////////////////////////////////////////////////////
    
    //open input and output link list files
    ifstream finput;
    finput.open(infile,fstream::in);
    
    int max_nodeID=0;
    int src, dest;
    int count = 0;;
    vector < int> nodeID_list;
    vector < vector <int> > links;
    
    
    while (!finput.eof()) {
        finput >> src >> dest ;
        if (max_nodeID < src | max_nodeID < dest ) {
            max_nodeID = max(src,dest);
            links.resize(max_nodeID+1);
        }
        links[src].push_back(dest);
        links[dest].push_back(src);
        count++;
    }
    count = count-1;
    
    finput.close();
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
    // run WTM model to generate WTM maps
    ////////////////////////////////////////////////////////////////////////////////
    
    vector < unsigned short> changed_nodes;// node ids for those that changed state
    vector < unsigned short> old_changed_nodes;// node ids for those that changed state
    
    vector < vector <int> > activation_times;// activation times
    vector < unsigned short> adopted;// boolean vector 1/0 of whos adopted
    vector < unsigned short> adopted_old;//vector < int> unadopt_node_list;// list of nodes to check = possible
    bool no_change;// determine if nothing has changed
    
    activation_times.resize(links.size());
    
    t = clock();
    
    //conduct filtrations starting at each node j
    for (int j=0; j<links.size(); ++j) {
        
        //set adopted to be all zeros
        adopted.resize(links.size());
        for (int i=0; i<links.size(); ++i)
            adopted[i] = 0;
        
        //set up activation times vector for seed at node j
        for (int i=0; i<links.size(); ++i)
            activation_times[j].push_back(links.size()*2);
        
        //make the cluster seed include node j and its neighbors
        activation_times[j][j] = 0;
        adopted[j] = 1;
        for (int neighb=0; neighb<links[j].size(); neighb++) {
            adopted[ links[j][neighb] ] = 1;
            activation_times[j][ links[j][neighb] ] = 0;
        }
        
        changed_nodes.clear();
        old_changed_nodes.clear();
        
        //add neighbors of contagion seed to boundary
        for (int neighb=0; neighb<links[j].size(); neighb++) {//node j's neighbors
            //  cerr << "    neighbor to check " << links[j][neighb] << endl;
            changed_nodes.push_back( links[j][neighb] );
        }
        
        //simulate WTM to get activation times
        for (int t=1; t<links.size(); ++t) {
            
            no_change = true;//assume there will be no change
            adopted_old = adopted;
            
            old_changed_nodes = changed_nodes;

            changed_nodes.clear();

            if (old_changed_nodes.size() > 0) {
                
                //go through nodes that changed stat during last timestep
                for (int k=0; k<old_changed_nodes.size(); ++k) {

                    //consider his neighbors
                    for (int l=0; l<links[ old_changed_nodes[k] ].size(); ++l) {
                        
                    int node2check = links[ old_changed_nodes[k] ][l];
                        
                        if (adopted[node2check]==0) {//if not adopted
                            
                            //computed number of neighbors that have adopted
                            int num_neighb_adopted = 0;
                            for (int neighb=0; neighb<links[ node2check ].size(); neighb++) {
                                num_neighb_adopted += adopted_old[ links[ node2check ][neighb]];
                            }
                            
                            //if fraction surpasses threshold
                            if ( ((double)num_neighb_adopted / (double)links[node2check].size() ) >T) {
                                adopted[node2check] = 1;//node adopts contagion
                                changed_nodes.push_back(node2check);//add node to set that changed this timestep
                                activation_times[j][node2check] = t;//save activation time
                                no_change = false;
                                
                            }
                            
                        }
                    }
                }
                
            }
            
        }
    }
    

    //save activation time matrix to a file
    ofstream foutput;
    foutput.open(outfile);
    for (int j=0; j<links.size(); ++j) {
        for (int i=0; i<links.size(); ++i) {
            foutput << activation_times[j][i];
            if (i<links.size()-1)
                foutput << "\t";
            else
                foutput << endl;
        }
    }
    foutput.close();

    
    
    return 0;
}

