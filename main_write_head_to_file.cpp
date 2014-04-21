#include "cal_tot_asa_of_res.h"

int main()
{
    
 
    ofstream file_distance("distance_of_four_atom.txt",ios::out);
    file_distance<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<"f_name"<<setw(12)<<"distance"<<endl;


    ofstream file0("count_of_total.txt",ios::out); 
    file0<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<"file/tot"<<setw(12)<<"total_c"<<endl;


    ofstream file1("count_of_target.txt",ios::out);

    int i;   
    string aa_list[20]={"ALA","VAL","PHE","PRO","MET","ILE","LEU","ASP","GLU","LYS",
                    "ARG","SER","THR","TYR","HIS","CYS","ASN","GLN","TRP","GLY"};
 
   
   file1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<"aa:";
   for(i=0;i<20;i++)
   {  
    file1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
   }
    file1<<endl;

   //consider neighbor 1 position
   ofstream file_n_1("count_of_n_1.txt",ios::out);
   file_n_1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<"aa:";
   for(i=0;i<20;i++)
   {  
    file_n_1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<aa_list[i];
   }
    file_n_1<<endl;

   //consider neighbor 2 position

   ofstream file_n_2("count_of_n_2.txt",ios::out);
   file_n_2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<"aa:";
   for(i=0;i<20;i++)
   {  
    file_n_2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<aa_list[i];
   }
    file_n_2<<endl;

   //consider neighbor 3 position

   ofstream file_n_3("count_of_n_3.txt",ios::out);
   file_n_3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(10)<<"aa:";
   for(i=0;i<20;i++)
   {  
    file_n_3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<aa_list[i];
   }
    file_n_3<<endl;


    ofstream file2("count_of_neighbor.txt",ios::out); 

    file2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(12)<<"file/nei"
        <<setw(6)<<"n_0"<<setw(6)<<"n_1"<<setw(6)<<"n_2"
        <<setw(6)<<"n_3"<<setw(6)<<"n_4"<<setw(6)<<"n_5"
        <<setw(6)<<"n_6"<<setw(6)<<"n_7"<<setw(6)<<"n_8"
        <<setw(6)<<"n_9"<<setw(6)<<"n_10"<<setw(6)<<"n_11"
        <<setw(6)<<"n_12"<<setw(6)<<"n_13"<<setw(6)<<"n_14"
        <<setw(6)<<"n_15"<<setw(6)<<"n_16"
        <<setw(6)<<"n_17"<<setw(6)<<"n_18"<<setw(6)<<"n_19"
        <<setw(6)<<"n_20"<<setw(6)<<"n_21"<<endl;

    ofstream file3("count_of_n_0_aa.txt",ios::out); 
    file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_0_aa:";
    for(i=0;i<20;i++)
    {  
    file3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file3<<endl;

    ofstream file4("count_of_n_1_aa.txt",ios::out); 
    file4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_1_aa:";
    for(i=0;i<20;i++)
    {  
    file4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file4<<endl;

    ofstream file5("count_of_n_2_aa.txt",ios::out); 
    file5<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_2_aa:";
    for(i=0;i<20;i++)
    {  
    file5<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file5<<endl;

    ofstream file6("count_of_n_3_aa.txt",ios::out); 
    file6<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_3_aa:";
    for(i=0;i<20;i++)
    {  
    file6<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file6<<endl;

    ofstream file7("count_of_n_4_aa.txt",ios::out); 
    file7<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_4_aa:";
    for(i=0;i<20;i++)
    {  
    file7<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file7<<endl;

    ofstream file8("count_of_n_5_aa.txt",ios::out); 
    file8<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_5_aa:";
    for(i=0;i<20;i++)
    {  
    file8<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file8<<endl;

    ofstream file9("count_of_n_6_aa.txt",ios::out); 
    file9<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_6_aa:";
    for(i=0;i<20;i++)
    {  
    file9<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file9<<endl;

    ofstream file10("count_of_n_7_aa.txt",ios::out); 
    file10<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_7_aa:";
    for(i=0;i<20;i++)
    {  
    file10<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file10<<endl;

    ofstream file11("count_of_n_8_aa.txt",ios::out); 
    file11<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_8_aa:";
    for(i=0;i<20;i++)
    {  
    file11<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file11<<endl;

    ofstream file12("count_of_n_9_aa.txt",ios::out); 
    file12<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(14)<<"file/n_9_aa:";
    for(i=0;i<20;i++)
    {  
    file12<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<aa_list[i];
    }
    file12<<endl;

    return 0;
}
