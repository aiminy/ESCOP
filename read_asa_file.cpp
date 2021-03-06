
#include "cal_tot_asa_of_res.h"
void read_asa_file(vector <aa> &w_c,coord_struct &asa, char *asa_file_name)
{
	ifstream asa_file;               // stream name for input *.asa file
	char dummy[RECORD_MAX],          // temporary C string for reading fields
	     tmp_id[FIELD_MAX];       	 // temporary C string, checks for "CA"
	register int  i;		 // fast access iterator

	asa.num_of_atoms=0;
	asa.num_of_res=0;
	
	asa_file.open(asa_file_name,ios::in);

	if (asa_file.bad()) 
	{
		cerr <<"Error: Unable to open "<<asa_file_name<<endl;
		exit (2);
	}//end of if

	while ((asa_file >> dummy)) 
	{
		if (strncmp(dummy,"ATOM",4)==0) 
		{
			++asa.num_of_atoms;
			asa_file >> dummy >> tmp_id;
			if (strncmp(tmp_id,"CA",2 )==0)
			++asa.num_of_res;
		}//end of if
		
		asa_file.getline(dummy,RECORD_MAX);
	}//end of while
	

	if (asa.num_of_atoms==0 || asa.num_of_res==0) 
	{
		cerr << "No usable residues in asa file.\n";
		exit (1);
	}// end of if

	asa.atom_num     = new int   [asa.num_of_atoms];
	asa.atom_id      = new char *[asa.num_of_atoms];
	asa.res_id       = new char *[asa.num_of_atoms];
	
	for (i=0; i<asa.num_of_atoms; ++i) 
	{
		asa.atom_id[i]   = new char [5];
		asa.res_id[i]    = new char [5];
	}
	
	asa.chain_id     = new char[asa.num_of_atoms+1];
	asa.res_num = new int    [asa.num_of_atoms];
	asa.x       = new double [asa.num_of_atoms];
	asa.y       = new double [asa.num_of_atoms];
	asa.z       = new double [asa.num_of_atoms];
	asa.asa     = new double [asa.num_of_atoms];

	asa_file.close();
	i=0;
	ifstream asa_file3;
        asa_file3.open(asa_file_name,ios::in);
	
	while (asa_file3>>dummy)
	{
		
		if (strncmp(dummy,"ATOM",4)!=0) 
			{
			asa_file3.getline(dummy,RECORD_MAX);
			}//end of if
		
		else {
                        asa_file3.ignore(2);
			asa_file3.get(dummy,6);
                        asa.atom_num[i]=atoi(dummy);
                        asa_file3.ignore(2);
			asa_file3.get(asa.atom_id[i],4);
                        asa_file3.ignore(1);
			asa_file3.get(asa.res_id[i],4);
                        asa_file3.ignore(1);
			asa_file3.get(asa.chain_id[i]);
			asa_file3.get(dummy,6); asa.res_num[i]=atoi(dummy);
			asa_file3.ignore(4);
			asa_file3.get(dummy,9); asa.x[i] = atof(dummy);
			asa_file3.get(dummy,9); asa.y[i] = atof(dummy);
			asa_file3.get(dummy,9); asa.z[i] = atof(dummy);
			asa_file3.get(dummy,9); 
                        //cout<<dummy<<endl;   
                        asa.asa[i]=atof(dummy);
			i++;
			asa_file3.getline(dummy,RECORD_MAX);
		     }//end of else
	}//end of while
	asa_file3.close();

        aa amino;
	amino.res_id       = new char [5];
        
	atom atom_of_res;
	atom_of_res.atom_id = new char [5];

	vector <aa> res; 
	int j;
	
	for (i=0;i<asa.num_of_atoms;i++)
	{
	  
	 if(res.empty())
	{
	   amino.chain_type=asa.chain_id[i];
	   amino.res_id=asa.res_id[i];
	   amino.res_num=asa.res_num[i];
	   res.push_back(amino);
	    
	for(j=0;j<asa.num_of_atoms;j++)
	{
	 if(*asa.res_id[j]==*res.back().res_id&&asa.res_num[j]==res.back().res_num&&asa.chain_id[j]==res.back().chain_type)
	   {
	   atom_of_res.atom_id=asa.atom_id[j];
	   atom_of_res.x=asa.x[j];
	   atom_of_res.y=asa.y[j];
	   atom_of_res.z=asa.z[j];
	   atom_of_res.asa=asa.asa[j];
	   amino.atom_str.push_back(atom_of_res);
	   }//end of if(*asa.re_id[j]..)
          	   
 	}//end of for(j=0)
	  w_c.push_back(amino);
         }// end of if (res.empty())

	 if(!res.empty())
	 {
          if(*asa.res_id[i]!=*res.back().res_id||asa.res_num[i]!=res.back().res_num||asa.chain_id[i]!=res.back().chain_type)
	  { 
	   amino.chain_type=asa.chain_id[i];
	   amino.res_id=asa.res_id[i];
	   amino.res_num=asa.res_num[i];
	   res.push_back(amino);
           
	   amino.atom_str.clear();
	   
	for(j=0;j<asa.num_of_atoms;j++)
	{
	 if(*asa.res_id[j]==*res.back().res_id&&asa.res_num[j]==res.back().res_num&&asa.chain_id[j]==res.back().chain_type)
	   {
	   atom_of_res.atom_id=asa.atom_id[j];
	   atom_of_res.x=asa.x[j];
	   atom_of_res.y=asa.y[j];
	   atom_of_res.z=asa.z[j];
	   atom_of_res.asa=asa.asa[j];
	   amino.atom_str.push_back(atom_of_res);
	   }//end of if(*asa.res_id[j]..)
          	   
 	}//end of for(j=0)
	  w_c.push_back(amino);
         }// end of if (*asa.res_id[i]!=)

	  } //end of if(!res.empty())

	}// end of for(i=0)

	return;

 }//end of function of read_asa_file
