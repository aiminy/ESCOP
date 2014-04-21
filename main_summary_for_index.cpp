#include "cal_tot_asa_of_res.h"

int main()
{
    
   int i,j,k,s;
   four_nu four_atom;
   vector<four_nu> std_nu;
   ifstream index_file;    

  for(s=0;s<20;s++)
  {
   for(i=0;i<9;i++)
   {
   for(j=0;j<9;j++)
    {
     for(k=0;k<9;k++)
     {
     four_atom.a=s;
     four_atom.b=i;
     four_atom.c=j;
     four_atom.d=k;
     std_nu.push_back(four_atom);
     }//end for
    }//end for 
   }//end for
  }//end fora

   for(i=0;i<std_nu.size();i++)
   {
    std_nu[i].p=0;
    std_nu[i].count=0;
    std_nu[i].angle=0; 
   }

   //cout<<"size"<<std_nu.size()<<endl;

   vector<four_nu> all_four_atom_chain;
   string t;
   int a,b,c,d; 
   double angle;
   index_file.open("index.txt",ios::in);

   while (index_file.get()!=-1)
   {
   index_file>>t;
   index_file>>a>>b>>c>>d;
   index_file>>angle;
   four_atom.a=a;
   four_atom.b=b;
   four_atom.c=c;
   four_atom.d=d;
   four_atom.angle=angle; 
   all_four_atom_chain.push_back(four_atom);
   index_file.get();
   }

    for(i=0;i<std_nu.size();i++)
    {
     cout<<std_nu[i].a<<" "
         <<std_nu[i].b<<" "
         <<std_nu[i].c<<" "
         <<std_nu[i].d<<endl;
    }


    double total;
    total=all_four_atom_chain.size();

   for(i=0;i<std_nu.size();i++)
   {
   for(j=0;j<all_four_atom_chain.size();j++)
    {
     if(all_four_atom_chain[j].a==std_nu[i].a&&
        all_four_atom_chain[j].b==std_nu[i].b&&
        all_four_atom_chain[j].c==std_nu[i].c&&
        all_four_atom_chain[j].d==std_nu[i].d)
      {
       std_nu[i].p=std_nu[i].p+1.0;
       std_nu[i].angle=std_nu[i].angle+all_four_atom_chain[j].angle;
       std_nu[i].count=std_nu[i].count+1;  
      }
    }
       if(std_nu[i].count!=0)
       {
       std_nu[i].angle=std_nu[i].angle/std_nu[i].count;
       }

       std_nu[i].p=std_nu[i].p/total;
   }

   ofstream outfile("index2.txt",ios::out);
   for(i=0;i<std_nu.size();i++)
   {
    if(std_nu[i].p>0)
    {
    outfile<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<"("
           <<setw(3)<<std_nu[i].a
           <<setw(3)<<std_nu[i].b
           <<setw(3)<<std_nu[i].c
           <<setw(3)<<std_nu[i].d<<setw(2)<<")"
           <<setiosflags(ios::fixed)
           <<setprecision(5)<<setw(12)<<std_nu[i].p
           <<setprecision(5)<<setw(12)<<std_nu[i].angle<<endl;
    }
   }


   //t_GLY_GLY_GLY
   ofstream outfile1("t_gly_gly_gly.txt",ios::out);
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].b==8&&all_four_atom_chain[j].c==8&&all_four_atom_chain[j].d==8)
      {
        outfile1<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<"("
                <<setw(3)<<all_four_atom_chain[j].a
                <<setw(3)<<all_four_atom_chain[j].b
                <<setw(3)<<all_four_atom_chain[j].c
                <<setw(3)<<all_four_atom_chain[j].d
                <<setprecision(5)<<setw(12)<<all_four_atom_chain[j].angle
                <<setw(2)<<")"
                <<endl;
      }
  }


    double total_angle1[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double count1[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double average1[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for(i=0;i<20;i++)
    {
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].a==i&&all_four_atom_chain[j].b==8&&all_four_atom_chain[j].c==8
           &&all_four_atom_chain[j].d==8)
      {
             total_angle1[i]=total_angle1[i]+all_four_atom_chain[j].angle;
             count1[i]=count1[i]+1;
      }
    }
     average1[i]=total_angle1[i]/count1[i];
   }
    
   ofstream outfile1_ave("t_gly_gly_gly_ave.txt",ios::out);

    for(i=0;i<20;i++)
     {   
     outfile1_ave<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<i
                 <<setprecision(5)<<setw(12)<<average1[i]<<endl;
    }


   //t_GLY_GLY_any
   ofstream outfile2("t_gly_gly_any.txt",ios::out);
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].b==8&&all_four_atom_chain[j].c==8)
      {
        outfile2<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<"("
                <<setw(3)<<all_four_atom_chain[j].a
                <<setw(3)<<all_four_atom_chain[j].b
                <<setw(3)<<all_four_atom_chain[j].c
                <<setw(3)<<all_four_atom_chain[j].d
                <<setprecision(5)<<setw(12)<<all_four_atom_chain[j].angle
                <<setw(2)<<")"
                <<endl;
      }
  }

    double total_angle2[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double count2[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double average2[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for(i=0;i<20;i++)
    {
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].a==i&&all_four_atom_chain[j].b==8&&all_four_atom_chain[j].c==8)
      {
             total_angle2[i]=total_angle2[i]+all_four_atom_chain[j].angle;
             count2[i]=count2[i]+1;
      }
    }
     average2[i]=total_angle2[i]/count2[i];
   }
    
   ofstream outfile2_ave("t_gly_gly_any_ave.txt",ios::out);

    for(i=0;i<20;i++)
     {   
     outfile2_ave<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<i
                 <<setprecision(5)<<setw(12)<<average2[i]<<endl;
    }

   //t_GLY_any_GLY
   ofstream outfile3("t_gly_any_gly.txt",ios::out);
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].b==8&&all_four_atom_chain[j].d==8)
      {
        outfile3<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<"("
                <<setw(3)<<all_four_atom_chain[j].a
                <<setw(3)<<all_four_atom_chain[j].b
                <<setw(3)<<all_four_atom_chain[j].c
                <<setw(3)<<all_four_atom_chain[j].d
                <<setprecision(5)<<setw(12)<<all_four_atom_chain[j].angle
                <<setw(2)<<")"
                <<endl;
      }
  }

    double total_angle3[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double count3[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double average3[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for(i=0;i<20;i++)
    {
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].a==i&&all_four_atom_chain[j].b==8&&all_four_atom_chain[j].d==8)
      {
             total_angle3[i]=total_angle3[i]+all_four_atom_chain[j].angle;
             count3[i]=count3[i]+1;
      }
    }
     average3[i]=total_angle3[i]/count3[i];
   }
    
   ofstream outfile3_ave("t_gly_any_gly_ave.txt",ios::out);

    for(i=0;i<20;i++)
     {   
     outfile3_ave<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<i
                 <<setprecision(5)<<setw(12)<<average3[i]<<endl;
    }

   //t_any_GLY_GLY
   ofstream outfile4("t_any_gly_gly.txt",ios::out);
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].c==8&&all_four_atom_chain[j].d==8)
      {
        outfile4<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<"("
                <<setw(3)<<all_four_atom_chain[j].a
                <<setw(3)<<all_four_atom_chain[j].b
                <<setw(3)<<all_four_atom_chain[j].c
                <<setw(3)<<all_four_atom_chain[j].d
                <<setprecision(5)<<setw(12)<<all_four_atom_chain[j].angle
                <<setw(2)<<")"
                <<endl;
      }
  }


    double total_angle4[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double count4[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double average4[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

        for(i=0;i<20;i++)
    {
   for(j=0;j<all_four_atom_chain.size();j++)
    {
        if(all_four_atom_chain[j].a==i&&all_four_atom_chain[j].c==8&&all_four_atom_chain[j].d==8)
      {
             total_angle4[i]=total_angle4[i]+all_four_atom_chain[j].angle;
             count4[i]=count4[i]+1;
      }
    }
     average4[i]=total_angle4[i]/count4[i];
   }
    
   ofstream outfile4_ave("t_any_gly_gly_ave.txt",ios::out);

    for(i=0;i<20;i++)
     {   
     outfile4_ave<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(2)<<i
                 <<setprecision(5)<<setw(12)<<average4[i]<<endl;
    }

    return 0;
}
