#include "cal_tot_asa_of_res.h"

void sum_for_closest_four_atom_chain(char *file_name,vector<target_many_neighbor> &closest_four_atom_chain)
{
        
    int i,j,t;
    char *temp;
    temp=new char[5];

    string aa_list[20]={"ALA","VAL","PHE","PRO","MET","ILE","LEU","ASP","GLU","LYS",
                    "ARG","SER","THR","TYR","HIS","CYS","ASN","GLN","TRP","GLY"};

      cout<<file_name<<endl;

     for(i=0;i<closest_four_atom_chain.size();i++)
        {
          cout<<closest_four_atom_chain[i].t.res_num<<" "<<closest_four_atom_chain[i].t.res_id
              <<" "<<closest_four_atom_chain[i].t.chain_type<<endl;
          for(j=0;j<3;j++)// 0: n1 1:n2 2:n3
          {
           //if(closest_four_atom_chain[i].n[j].d!=0)
           //{
           cout<<closest_four_atom_chain[i].n[j].res_num<<" "<<closest_four_atom_chain[i].n[j].res_id<<" "
               <<closest_four_atom_chain[i].n[j].chain_type<<" "
               <<closest_four_atom_chain[i].n[j].d<<endl;
           //}  
          }
        }

        ofstream file_distance("distance_of_four_atom.txt",ios::out|ios::app);
	for(i=0;i<closest_four_atom_chain.size();i++)
        {
          for(j=0;j<3;j++)
          {
           //if(closest_four_atom_chain[i].n[j].d>3)
           //{
            file_distance<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<file_name
                         <<setw(12)<<closest_four_atom_chain[i].n[j].d<<endl;
           //}  
          }
        }


         double  average_n_d;   
        for(i=0;i<closest_four_atom_chain.size();i++)
        {
        //  cout<<t_m_n_chain[i].t.res_num<<" "<<t_m_n_chain[i].t.res_id<<endl;
          average_n_d=0;  
          for(j=0;j<3;j++)
          {
           average_n_d=average_n_d+closest_four_atom_chain[i].n[j].d;  
          }
          closest_four_atom_chain[i].t.ave_n_d=average_n_d/3.0;
        }

        ofstream file_average_distance("average_distance_of_neighbors.txt",ios::out|ios::app);
	for(i=0;i<closest_four_atom_chain.size();i++)
        {
            file_average_distance<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<file_name
                         <<setw(12)<<closest_four_atom_chain[i].t.ave_n_d<<endl;
        }

   double total_count;
   total_count=closest_four_atom_chain.size();

   ofstream file0("count_of_total.txt",ios::out|ios::app); 
   file0<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<file_name<<setw(12)<<total_count<<endl;

   double n_t[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].t.res_id,"ALA")==0)
     {n_t[0]=n_t[0]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"VAL")==0)
     {n_t[1]=n_t[1]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"PHE")==0)
     {n_t[2]=n_t[2]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"PRO")==0)
     {n_t[3]=n_t[3]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"MET")==0)
     {n_t[4]=n_t[4]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"ILE")==0)
     {n_t[5]=n_t[5]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"LEU")==0)
     {n_t[6]=n_t[6]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"ASP")==0)
     {n_t[7]=n_t[7]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"GLU")==0)
     {n_t[8]=n_t[8]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"LYS")==0)
     {n_t[9]=n_t[9]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"ARG")==0)
     {n_t[10]=n_t[10]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"SER")==0)
     {n_t[11]=n_t[11]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"THR")==0)
     {n_t[12]=n_t[12]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"TYR")==0)
     {n_t[13]=n_t[13]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"HIS")==0)
     {n_t[14]=n_t[14]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"CYS")==0)
     {n_t[15]=n_t[15]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"ASN")==0)
     {n_t[16]=n_t[16]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"GLN")==0)
     {n_t[17]=n_t[17]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"TRP")==0)
     {n_t[18]=n_t[18]+1;}
   else if(strcmp(closest_four_atom_chain[i].t.res_id,"GLY")==0)
     {n_t[19]=n_t[19]+1;}
   else 
     {n_t[20]=n_t[20]+1;}
}


   ofstream file1("count_of_target.txt",ios::out|ios::app); 
   file1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<file_name;
   
   for(i=0;i<20;i++)
   {  
   file1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<n_t[i];
   }
    file1<<endl;

   //consider neighbor 1 position
   double n_1[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ALA")==0)
     {n_1[0]=n_1[0]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"VAL")==0)
     {n_1[1]=n_1[1]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"PHE")==0)
     {n_1[2]=n_1[2]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"PRO")==0)
     {n_1[3]=n_1[3]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"MET")==0)
     {n_1[4]=n_1[4]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ILE")==0)
     {n_1[5]=n_1[5]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"LEU")==0)
     {n_1[6]=n_1[6]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ASP")==0)
     {n_1[7]=n_1[7]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"GLU")==0)
     {n_1[8]=n_1[8]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0)
     {n_1[9]=n_1[9]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0)
     {n_1[10]=n_1[10]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"SER")==0)
     {n_1[11]=n_1[11]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"THR")==0)
     {n_1[12]=n_1[12]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"TYR")==0)
     {n_1[13]=n_1[13]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)
     {n_1[14]=n_1[14]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"CYS")==0)
     {n_1[15]=n_1[15]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ASN")==0)
     {n_1[16]=n_1[16]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"GLN")==0)
     {n_1[17]=n_1[17]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"TRP")==0)
     {n_1[18]=n_1[18]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"GLY")==0)
     {n_1[19]=n_1[19]+1;}
   else 
     {n_1[20]=n_1[20]+1;}
}


   ofstream file_n_1("count_of_n_1.txt",ios::out|ios::app); 
   file_n_1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<file_name;
   
   for(i=0;i<20;i++)
   {  
   file_n_1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<n_1[i];
   }
    file_n_1<<endl;


    //consider neighbor 2 position

   double n_2[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ALA")==0)
     {n_2[0]=n_2[0]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"VAL")==0)
     {n_2[1]=n_2[1]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"PHE")==0)
     {n_2[2]=n_2[2]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"PRO")==0)
     {n_2[3]=n_2[3]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"MET")==0)
     {n_2[4]=n_2[4]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ILE")==0)
     {n_2[5]=n_2[5]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"LEU")==0)
     {n_2[6]=n_2[6]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0)
     {n_2[7]=n_2[7]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)
     {n_2[8]=n_2[8]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"LYS")==0)
     {n_2[9]=n_2[9]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ARG")==0)
     {n_2[10]=n_2[10]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"SER")==0)
     {n_2[11]=n_2[11]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"THR")==0)
     {n_2[12]=n_2[12]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"TYR")==0)
     {n_2[13]=n_2[13]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"HIS")==0)
     {n_2[14]=n_2[14]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"CYS")==0)
     {n_2[15]=n_2[15]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0)
     {n_2[16]=n_2[16]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)
     {n_2[17]=n_2[17]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"TRP")==0)
     {n_2[18]=n_2[18]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"GLY")==0)
     {n_2[19]=n_2[19]+1;}
   else 
     {n_2[20]=n_2[20]+1;}
}


   ofstream file_n_2("count_of_n_2.txt",ios::out|ios::app); 
   file_n_2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<file_name;
   
   for(i=0;i<20;i++)
   {  
   file_n_2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<n_2[i];
   }
    file_n_2<<endl;


   //consider neighbor 3 position
   double n_3[21]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ALA")==0)
     {n_3[0]=n_3[0]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"VAL")==0)
     {n_3[1]=n_3[1]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"PHE")==0)
     {n_3[2]=n_3[2]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"PRO")==0)
     {n_3[3]=n_3[3]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"MET")==0)
     {n_3[4]=n_3[4]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ILE")==0)
     {n_3[5]=n_3[5]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"LEU")==0)
     {n_3[6]=n_3[6]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ASP")==0)
     {n_3[7]=n_3[7]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"GLU")==0)
     {n_3[8]=n_3[8]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"LYS")==0)
     {n_3[9]=n_3[9]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ARG")==0)
     {n_3[10]=n_3[10]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0)
     {n_3[11]=n_3[11]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0)
     {n_3[12]=n_3[12]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"TYR")==0)
     {n_3[13]=n_3[13]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"HIS")==0)
     {n_3[14]=n_3[14]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"CYS")==0)
     {n_3[15]=n_3[15]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ASN")==0)
     {n_3[16]=n_3[16]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"GLN")==0)
     {n_3[17]=n_3[17]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"TRP")==0)
     {n_3[18]=n_3[18]+1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"GLY")==0)
     {n_3[19]=n_3[19]+1;}
   else 
     {n_3[20]=n_3[20]+1;}
}

   ofstream file_n_3("count_of_n_3.txt",ios::out|ios::app); 
   file_n_3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<file_name;
   
   for(i=0;i<20;i++)
   {  
   file_n_3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<n_3[i];
   }
    file_n_3<<endl;



 for(j=0;j<20;j++)
 { 
   strcpy(temp,aa_list[j].c_str());
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].t.res_id,temp)==0)
   {
       closest_four_atom_chain[i].t_id=j;
   }//end if
   }//end for 
}//end for

   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)
     {closest_four_atom_chain[i].n1_id=0;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ASP")==0||
           strcmp(closest_four_atom_chain[i].n[1].res_id,"GLU")==0)
     {closest_four_atom_chain[i].n1_id=1;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"GLN")==0)
     {closest_four_atom_chain[i].n1_id=2;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"SER")==0||
           strcmp(closest_four_atom_chain[i].n[1].res_id,"THR")==0)
     {closest_four_atom_chain[i].n1_id=3;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"ALA")==0)
     {closest_four_atom_chain[i].n1_id=4;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"CYS")==0)
     {closest_four_atom_chain[i].n1_id=5;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"PRO")==0)
     {closest_four_atom_chain[i].n1_id=6;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"VAL")==0)
     {closest_four_atom_chain[i].n1_id=7;}
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"GLY")==0)
     {closest_four_atom_chain[i].n1_id=8;}
   }//for

   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"HIS")==0)
     {closest_four_atom_chain[i].n2_id=0;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
           strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)
     {closest_four_atom_chain[i].n2_id=1;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)
     {closest_four_atom_chain[i].n2_id=2;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"SER")==0||
           strcmp(closest_four_atom_chain[i].n[2].res_id,"THR")==0)
     {closest_four_atom_chain[i].n2_id=3;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"ALA")==0)
     {closest_four_atom_chain[i].n2_id=4;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"CYS")==0)
     {closest_four_atom_chain[i].n2_id=5;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"PRO")==0)
     {closest_four_atom_chain[i].n2_id=6;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"VAL")==0)
     {closest_four_atom_chain[i].n2_id=7;}
   else if(strcmp(closest_four_atom_chain[i].n[2].res_id,"GLY")==0)
     {closest_four_atom_chain[i].n2_id=8;}
   }//for

   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"HIS")==0)
     {closest_four_atom_chain[i].n3_id=0;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ASP")==0||
           strcmp(closest_four_atom_chain[i].n[3].res_id,"GLU")==0)
     {closest_four_atom_chain[i].n3_id=1;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLN")==0)
     {closest_four_atom_chain[i].n3_id=2;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
           strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0)
     {closest_four_atom_chain[i].n3_id=3;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"ALA")==0)
     {closest_four_atom_chain[i].n3_id=4;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"CYS")==0)
     {closest_four_atom_chain[i].n3_id=5;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"PRO")==0)
     {closest_four_atom_chain[i].n3_id=6;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"VAL")==0)
     {closest_four_atom_chain[i].n3_id=7;}
   else if(strcmp(closest_four_atom_chain[i].n[3].res_id,"GLY")==0)
     {closest_four_atom_chain[i].n3_id=8;}
   }//for


   four_nu four_atom;
   vector<four_nu> four_atom_chain;
   for(i=0;i<closest_four_atom_chain.size();i++)
   {
    four_atom.a=closest_four_atom_chain[i].t_id;
    four_atom.b=closest_four_atom_chain[i].n1_id;
    four_atom.c=closest_four_atom_chain[i].n2_id;
    four_atom.d=closest_four_atom_chain[i].n3_id;
    four_atom.angle=closest_four_atom_chain[i].t.angle;
    four_atom_chain.push_back(four_atom);    
   }

    ofstream file_index("index.txt",ios::out|ios::app); 
   
   for(i=0;i<four_atom_chain.size();i++)
   {
    file_index<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<file_name
              <<setw(6)<<four_atom_chain[i].a
              <<setw(6)<<four_atom_chain[i].b
              <<setw(6)<<four_atom_chain[i].c
              <<setw(6)<<four_atom_chain[i].d
              <<setw(6)<<four_atom_chain[i].angle<<endl;
   }    

 //consider to count the combination of target residue and its neighbors  
 double c[20][10];

 for(i=0;i<20;i++)
 {
  for(j=0;j<10;j++)
  {
   c[i][j]=0;
  }
 }
 

 for(j=0;j<20;j++)
 { 

   strcpy(temp,aa_list[j].c_str());
   for(i=0;i<closest_four_atom_chain.size();i++)
   {

   if(strcmp(closest_four_atom_chain[i].t.res_id,temp)==0)
   {
   if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
       strcmp(closest_four_atom_chain[i].n[2].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"HIS")==0)&&(
       strcmp(closest_four_atom_chain[i].n[3].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"HIS")==0))
        {
        c[j][0]=c[j][0]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLU")==0))
        {
        c[j][1]=c[j][1]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLN")==0))
        {
        c[j][2]=c[j][2]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"THR")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"THR")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {
        c[j][3]=c[j][3]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ALA")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ALA")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ALA")==0))
        {
       //   cout<<"n1"<<closest_four_atom_chain[i].n[1].res_id<<endl;
       //   cout<<"n1"<<closest_four_atom_chain[i].n[2].res_id<<endl;
       //   cout<<"n1"<<closest_four_atom_chain[i].n[3].res_id<<endl;
        c[j][4]=c[j][4]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"CYS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"CYS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"CYS")==0))
        {
        c[j][5]=c[j][5]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"PRO")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"PRO")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PRO")==0))
        {
        c[j][6]=c[j][6]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"VAL")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"VAL")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"VAL")==0))
        {
         c[j][7]=c[j][7]+1; 
        } //end if
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"GLY")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLY")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLY")==0))
        {
        c[j][8]=c[j][8]+1; 
        }
    else
        {c[j][9]=c[j][9]+1;}

    } //end first if
}//for
}//for

/*
     cout<<num_0[0]<<endl;
     cout<<num_1[0]<<endl;
     cout<<num_2[0]<<endl;
     cout<<num_3[0]<<endl;
     cout<<num_4[0]<<endl;
     cout<<num_5[0]<<endl;
     cout<<num_6[0]<<endl;
     cout<<num_7[0]<<endl;
     cout<<num_8[0]<<endl;
     cout<<num_9[0]<<endl;
*/
     
   ofstream file_n_0_aa("count_of_n_0_aa.txt",ios::out|ios::app); 
   file_n_0_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   
   for(i=0;i<20;i++)
   {  
    file_n_0_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][0];
   }
    file_n_0_aa<<endl;

   ofstream file_n_1_aa("count_of_n_1_aa.txt",ios::out|ios::app); 
   file_n_1_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;

   for(i=0;i<20;i++)
   {  
    file_n_1_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][1];
   }
    file_n_1_aa<<endl;

   ofstream file_n_2_aa("count_of_n_2_aa.txt",ios::out|ios::app);
   file_n_2_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;

   for(i=0;i<20;i++)
   {  
    file_n_2_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][2];
   }
    file_n_2_aa<<endl;

   ofstream file_n_3_aa("count_of_n_3_aa.txt",ios::out|ios::app); 
   file_n_3_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_3_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][3];
   }
    file_n_3_aa<<endl;

   ofstream file_n_4_aa("count_of_n_4_aa.txt",ios::out|ios::app); 
   file_n_4_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;

   for(i=0;i<20;i++)
   {  
    file_n_4_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][4];
   }
    file_n_4_aa<<endl;

   ofstream file_n_5_aa("count_of_n_5_aa.txt",ios::out|ios::app);
   file_n_5_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_5_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][5];
   }
    file_n_5_aa<<endl;

   ofstream file_n_6_aa("count_of_n_6_aa.txt",ios::out|ios::app);
   file_n_6_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_6_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][6];
   }
    file_n_6_aa<<endl;

   ofstream file_n_7_aa("count_of_n_7_aa.txt",ios::out|ios::app);
   file_n_7_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_7_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][7];
   }
    file_n_7_aa<<endl;

   ofstream file_n_8_aa("count_of_n_8_aa.txt",ios::out|ios::app);
   file_n_8_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_8_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][8];
   }
    file_n_8_aa<<endl;

   ofstream file_n_9_aa("count_of_n_9_aa.txt",ios::out|ios::app);
   file_n_9_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<file_name;
   for(i=0;i<20;i++)
   {  
    file_n_9_aa<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<c[i][9];
   }
    file_n_9_aa<<endl;



   //consider to count three neighbors only  
   double n[22]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   for(i=0;i<closest_four_atom_chain.size();i++)
   {
   if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
       strcmp(closest_four_atom_chain[i].n[2].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[2].res_id,"HIS")==0)&&(
       strcmp(closest_four_atom_chain[i].n[3].res_id,"ARG")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"LYS")==0||
       strcmp(closest_four_atom_chain[i].n[3].res_id,"HIS")==0))
        {
        n[0]=n[0]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLU")==0))
        {
        n[1]=n[1]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLN")==0))
        {
        n[2]=n[2]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"THR")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"THR")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {
        n[3]=n[3]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ALA")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ALA")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ALA")==0))
        {
        n[4]=n[4]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"CYS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"CYS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"CYS")==0))
        {
        n[5]=n[5]+1; 
        }
   else if(strcmp(closest_four_atom_chain[i].n[1].res_id,"PRO")==0&&
           strcmp(closest_four_atom_chain[i].n[2].res_id,"PRO")==0&&
           strcmp(closest_four_atom_chain[i].n[3].res_id,"PRO")==0)
        {
        n[6]=n[6]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"VAL")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"VAL")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"VAL")==0))
        {
         n[7]=n[7]+1; 
        } //end if
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"GLY")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLY")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLY")==0))
        {
        n[8]=n[8]+1; 
        }
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_2
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLN")==0))
        {n[9]=n[9]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_3
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[10]=n[10]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_4
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ALA")==0))
        {n[11]=n[11]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_5
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"CYS")==0))
        {n[12]=n[12]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_6
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PRO")==0))
        {n[13]=n[13]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_7
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TRP")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"TYR")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"PHE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"MET")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"LEU")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"ILE")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"VAL")==0))
        {n[14]=n[14]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_1_8
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASP")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLU")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"GLY")==0))
        {n[15]=n[15]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_3
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[16]=n[16]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_4
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[17]=n[17]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_5
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[18]=n[18]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_6
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[19]=n[19]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_7
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[20]=n[20]+1;}
   else if((strcmp(closest_four_atom_chain[i].n[1].res_id,"ARG")==0|| //0_2_8
            strcmp(closest_four_atom_chain[i].n[1].res_id,"LYS")==0||
            strcmp(closest_four_atom_chain[i].n[1].res_id,"HIS")==0)&&(
            strcmp(closest_four_atom_chain[i].n[2].res_id,"ASN")==0||
            strcmp(closest_four_atom_chain[i].n[2].res_id,"GLN")==0)&&(
            strcmp(closest_four_atom_chain[i].n[3].res_id,"SER")==0||
            strcmp(closest_four_atom_chain[i].n[3].res_id,"THR")==0))
        {n[21]=n[21]+1;}
          
}//for

     
    ofstream file2("count_of_neighbor.txt",ios::out|ios::app); 

    file2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<file_name
        <<setw(6)<<n[0]<<setw(6)<<n[1]<<setw(6)<<n[2]
        <<setw(6)<<n[3]<<setw(6)<<n[4]<<setw(6)<<n[5]
        <<setw(6)<<n[6]<<setw(6)<<n[7]<<setw(6)<<n[8]
        <<setw(6)<<n[9]<<setw(6)<<n[10]<<setw(6)<<n[11]
        <<setw(6)<<n[12]<<setw(6)<<n[13]<<setw(6)<<n[14]
        <<setw(6)<<n[15]<<setw(6)<<n[16]
        <<setw(6)<<n[17]<<setw(6)<<n[18]<<setw(6)<<n[19]
        <<setw(6)<<n[20]<<setw(6)<<n[21]<<endl;

/*
   ofstream file3("summary_of_count.txt",ios::out|ios::app);
  
   file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"t/n"
        <<setw(6)<<"c_0"<<setw(6)<<"c_1"<<setw(6)<<"c_2"
        <<setw(6)<<"c_3"<<setw(6)<<"c_4"<<setw(6)<<"c_5"
        <<setw(6)<<"c_6"<<setw(6)<<"c_7"<<setw(6)<<"c_8"
        <<setw(6)<<"c_9"<<setw(6)<<"total"<<endl;
   
   for(i=0;i<20;i++)
   {  
   file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i]
        <<setw(6)<<c[i][0]<<setw(6)<<c[i][1]<<setw(6)<<c[i][2]
        <<setw(6)<<c[i][3]<<setw(6)<<c[i][4]<<setw(6)<<c[i][5]
        <<setw(6)<<c[i][6]<<setw(6)<<c[i][7]<<setw(6)<<c[i][8]
        <<setw(6)<<c[i][9]<<setw(6)<<n_t[i]<<endl;
   }

//   file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"other"
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<n_t[20]<<endl;
   
   file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"total"
        <<setw(6)<<n[0]<<setw(6)<<n[1]<<setw(6)<<n[2]
        <<setw(6)<<n[3]<<setw(6)<<n[4]<<setw(6)<<n[5]
        <<setw(6)<<n[6]<<setw(6)<<n[7]<<setw(6)<<n[8]
        <<setw(6)<<n[9]<<setw(6)<<total_count<<endl;
  
   double p[20][10];
  
   for(i=0;i<20;i++)
    {
      for(j=0;j<10;j++)
      {p[i][j]=0;}
    }

   for(i=0;i<20;i++)
    {
      for(j=0;j<10;j++)

      {if(n[j]!=0&&n_t[i]!=0)
       {
       p[i][j]=(c[i][j])/(n[j]*n_t[i]);
       }
       else
       {p[i][j]=0;}
      }
    }


   ofstream file4("summary_of_probability.txt",ios::out|ios::app);
  
   file4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"t/n"
        <<setw(12)<<"c_0"<<setw(12)<<"c_1"<<setw(12)<<"c_2"
        <<setw(12)<<"c_3"<<setw(12)<<"c_4"<<setw(12)<<"c_5"
        <<setw(12)<<"c_6"<<setw(12)<<"c_7"<<setw(12)<<"c_8"
        <<setw(12)<<"c_9"<<setw(12)<<"total"<<endl;
   
   for(i=0;i<20;i++)
   {  
   file4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i]
        <<setiosflags(ios::scientific)
        <<setprecision(5)<<setw(12)<<p[i][0]
        <<setprecision(5)<<setw(12)<<p[i][1]
        <<setprecision(5)<<setw(12)<<p[i][2]
        <<setprecision(5)<<setw(12)<<p[i][3]
        <<setprecision(5)<<setw(12)<<p[i][4]
        <<setprecision(5)<<setw(12)<<p[i][5]
        <<setprecision(5)<<setw(12)<<p[i][6]
        <<setprecision(5)<<setw(12)<<p[i][7]
        <<setprecision(5)<<setw(12)<<p[i][8]
        <<setprecision(5)<<setw(12)<<p[i][9]
        <<setprecision(5)<<setw(12)<<n_t[i]<<endl;
   }

//   file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"other"
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<""<<setw(6)<<""
//        <<setw(6)<<""<<setw(6)<<n_t[20]<<endl;
   
   file4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"total"
        <<setw(12)<<n[0]<<setw(12)<<n[1]<<setw(12)<<n[2]
        <<setw(12)<<n[3]<<setw(12)<<n[4]<<setw(12)<<n[5]
        <<setw(12)<<n[6]<<setw(12)<<n[7]<<setw(12)<<n[8]
        <<setw(12)<<n[9]<<setw(12)<<total_count<<endl;
*/

}//end of sum
