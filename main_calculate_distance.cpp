//@Aimin Yan
//filename:read_asa.cpp
//usage: ./read_asa *.asa
//read monomer protein file *.asa in this directory calculated from naccess
//calculate omega angle and indentify SASA of residue and atom
//output to file angle_asa.txt
//angle_asa.txt format:
//res_num res_name chain_id res_omega res_rel_asa
//for example:
//  1       GLU       A       67.00     34.87 
 
#include "cal_tot_asa_of_res.h"

int main (int argc, char *argv[]) {
	
	coord_struct  asa;
	int    i,j;
	char p_n[20];
	vector<char> chain;
        vector<aa> w_c,sp,w_c2;	
        vector<aa_asa> std_asa; 
        vector<naa2>w_c4;
        vector<naa>w_c3;

        vector<residue>a_w_c;
	
	//vector<complex_aa> total_aa;   
	
	char *c_a;
	c_a=new char[1];
	

         vector<protein> whole_protein;
         vector<aa> wc_c_aa2,wc_t_aa2,apc;
         vector<naa> wc_c_aa3,wc_t_aa3;
         //vector<atom> four_body;
         //vector< vector<atom> > closest_four_atom_chain; 
         vector<target_many_neighbor> t_m_n_chain_for_all_protein;
	
        if(argc<2)
        {cout<<"usage:"<<"./main_calculate_distance asa_file"<<endl;}

	 //ifstream protein_file(argv[1],ios::in);

         //char *file_name;
         //file_name=new char[10];

         //while(protein_file>>file_name)
        //{
        //cout<<file_name<<endl;
	read_asa_file(w_c,asa,argv[1]);
        //calculate_three_nearest_neighbors(w_c,sp);
        inter_atom_distance_matrix(w_c,sp,t_m_n_chain_for_all_protein);
         
        //regenerate_pdb_file(w_c);   
        //output_global_angle(w_c);

        //sp.clear(); 
        //w_c.clear();  
        //}//end of while 
         //end of read different protein 

          //cout<<"total triangle:"<<closest_four_atom_chain.size()<<endl;

         //protein_file.close(); 
            
        sum_for_closest_four_atom_chain(argv[1],t_m_n_chain_for_all_protein);       
        // sum_for_local(sp);

        return(0);
     }//end of main
