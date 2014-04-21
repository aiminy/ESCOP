#include "cal_tot_asa_of_res.h"
#include "tmatrix.h"
#include "Vector3D.h"

template <class T>

//int partition(vector<double> &, int, int);
//void quicksort(vector<double> &, int, int);

int find_minimum(vector<neighbor> &,int,int);
void swap(neighbor &a, neighbor &b);
void sort(vector<neighbor> &,int);

double angle_between_two_vectors(Vector3D A,Vector3D B);

double angle_between_two_vectors(Vector3D A,Vector3D B)
{
  double A_DOT_B=A.DotProduct(B);
  double A_length=A.CalculateLength();
  double B_length=B.CalculateLength();
  double AB=A_length*B_length;
  double angle=180*acos(A_DOT_B/AB)/3.14;
  return angle;
}

/*
int partition(vector <atom_pair> &array, int top, int bottom)
{
     atom_pair x1;

     x1=array[top];    
     int i = top - 1;
     int j = bottom + 1;
     atom_pair temp;

     do
     {
           do      
           {
           j--;
           }while (x1.d>array[j].d);

          do  
          {
          i++;
          } while (x1.d <array[i].d);

          if (i < j)
         { 
          temp = array[i];    // switch elements at positions i and j
          //assign_atom_pair(temp,array[i]); 
          array[i] = array[j];
          //assign_atom_pair(array[i],array[j]); 
          array[j] = temp;
          //assign_atom_pair(array[j],temp); 
         }

     }while (i < j);
     
     return j;           // returns middle index 
}

void quicksort(vector <atom_pair> &array, int top, int bottom)
{
      // top = subscript of beginning of vector being considered
      // bottom = subscript of end of vector being considered
      // this process uses recursion - the process of calling itself
     int middle;
     if (top < bottom)
    {
          middle = partition(array, top, bottom);
          quicksort(array, top, middle);   // sort top partition
          quicksort(array, middle+1, bottom);    // sort bottom partition
     }
     return;
}
*/

int find_minimum(vector<neighbor> &n_l,int start,int end)
{
 int i,location;
 neighbor minimum;
 minimum.d=n_l[start].d;
 minimum.res_num=n_l[start].res_num;
 minimum.res_id=new char[5];
 minimum.res_id=n_l[start].res_id;
 minimum.chain_type=n_l[start].chain_type;

 
 location=start;
 
 for(i=start+1;i<=end;i++)
 {
  if(n_l[i].d<minimum.d)
   {
     minimum.d=n_l[i].d;
     minimum.res_num=n_l[i].res_num;
     minimum.res_id=new char[5];
     minimum.res_id=n_l[i].res_id;
     minimum.chain_type=n_l[i].chain_type;
     location=i;
   }
 } 
 return location;
}


void swap(neighbor &a, neighbor &b)
{
  neighbor temp;
  temp.res_id=new char[5];
 
  temp.d=a.d;
  temp.res_num=a.res_num;
  temp.res_id=a.res_id;
  temp.chain_type=a.chain_type;



  a.d=b.d;
  a.res_num=b.res_num;
  a.res_id=new char[5];
  a.res_id=b.res_id;
  a.chain_type=b.chain_type;


  b.d=temp.d;
  b.res_num=temp.res_num;
  b.res_id=new char[5];
  b.res_id=temp.res_id;
  b.chain_type=temp.chain_type;


}
void sort(vector<neighbor> &array,int size)
{
   int i,location;
   for(i=0;i<=array.size()-2;i++)
  {
    location=find_minimum(array,i,array.size()-1);
    swap(array[i],array[location]);
  }
}


void inter_atom_distance_matrix(vector<aa> &w_c,vector<aa> &sp,vector<target_many_neighbor> &t_m_n_chain)
{
        int i,j,rows,cols,num_total_atom,num_of_calpha,rows1,cols1;
            

         num_total_atom=0;
	
	for(i=0;i<w_c.size();i++)
	 {
          num_total_atom=num_total_atom+w_c[i].atom_str.size();
         }

         rows=cols=num_total_atom;
  
       
        vector <atom> chain_all_atom;
 
        for(i=0;i<w_c.size();i++)
         {
          for(j=0;j<w_c[i].atom_str.size();j++)
           {
            w_c[i].atom_str[j].res_num=w_c[i].res_num;
            w_c[i].atom_str[j].chain_type=w_c[i].chain_type;
            w_c[i].atom_str[j].res_id=new char[5];
            w_c[i].atom_str[j].res_id=w_c[i].res_id;
            chain_all_atom.push_back(w_c[i].atom_str[j]);
           }
         }

          /*
          for(i=0;i<chain_all_atom.size();i++)
          { 
          if(strcmp(chain_all_atom[i].atom_id,"CA ")==0)
           {cout<<chain_all_atom[i].res_num<<" "<<chain_all_atom[i].res_id<<" "<<chain_all_atom[i].atom_id
                <<chain_all_atom[i].asa<<" "<<chain_all_atom[i].chain_type<<endl;}
          }
           */
 


        //center of structure
/*
          double ws_x,ws_y,ws_z;
          double how_many_atom;      
          ws_x=0; 
          ws_y=0; 
          ws_z=0; 
          how_many_atom=0;

       for(i=0;i<w_c.size();i++)
         {
                  w_c[i].global_angle=0;  
            
          for(j=0;j<w_c[i].atom_str.size();j++)
           {
            ws_x=ws_x+w_c[i].atom_str[j].x;
            ws_y=ws_y+w_c[i].atom_str[j].y;
            ws_z=ws_z+w_c[i].atom_str[j].z;
            how_many_atom=how_many_atom+1;
           }
         }

            ws_x=ws_x/how_many_atom;
            ws_y=ws_y/how_many_atom;
            ws_z=ws_z/how_many_atom;

         // cout<<ws_x<<" "<<ws_y<<" "<<ws_z<<" "<<endl;   
*/

         vector<atom> surface_calpha_atom;

          //include GLY
         for(i=0;i<chain_all_atom.size();i++)
          {
           if(strcmp(chain_all_atom[i].atom_id,"CA ")==0&&chain_all_atom[i].asa<=5.0000&&
               (
                strncmp(chain_all_atom[i].res_id,"ALA",3)==0||
                strncmp(chain_all_atom[i].res_id,"VAL",3)==0||
                strncmp(chain_all_atom[i].res_id,"PHE",3)==0||
                strncmp(chain_all_atom[i].res_id,"PRO",3)==0||
                strncmp(chain_all_atom[i].res_id,"MET",3)==0||
                strncmp(chain_all_atom[i].res_id,"ILE",3)==0||
                strncmp(chain_all_atom[i].res_id,"LEU",3)==0||
                strncmp(chain_all_atom[i].res_id,"ASP",3)==0||
                strncmp(chain_all_atom[i].res_id,"GLU",3)==0||
                strncmp(chain_all_atom[i].res_id,"LYS",3)==0||
                strncmp(chain_all_atom[i].res_id,"ARG",3)==0||
                strncmp(chain_all_atom[i].res_id,"SER",3)==0||
                strncmp(chain_all_atom[i].res_id,"THR",3)==0||
                strncmp(chain_all_atom[i].res_id,"TYR",3)==0||
                strncmp(chain_all_atom[i].res_id,"HIS",3)==0||
                strncmp(chain_all_atom[i].res_id,"CYS",3)==0||
                strncmp(chain_all_atom[i].res_id,"ASN",3)==0||
                strncmp(chain_all_atom[i].res_id,"GLN",3)==0||
                strncmp(chain_all_atom[i].res_id,"TRP",3)==0||
                strncmp(chain_all_atom[i].res_id,"GLY",3)==0
               ))
             {
              //chain_all_atom[i].atomID=chain_all_atom[i].res_num; 
              surface_calpha_atom.push_back(chain_all_atom[i]);
              //k++;
             }//end if
          }//end for

         
         double d;

         target_many_neighbor t_m_n;
         neighbor new_n; 
         //vector<target_many_neighbor> t_m_n_chain;
 
         for(i=0;i<surface_calpha_atom.size();i++)
         {

            t_m_n.t.res_num=surface_calpha_atom[i].res_num;
            t_m_n.t.res_id=new char[5];
            t_m_n.t.res_id=surface_calpha_atom[i].res_id;
            t_m_n.t.chain_type=surface_calpha_atom[i].chain_type;
            t_m_n.t.x=surface_calpha_atom[i].x;
            t_m_n.t.y=surface_calpha_atom[i].y;
            t_m_n.t.z=surface_calpha_atom[i].z;


           for(j=0;j<surface_calpha_atom.size();j++)
            {
            if(surface_calpha_atom[j].res_num==t_m_n.t.res_num&&
               surface_calpha_atom[j].chain_type==t_m_n.t.chain_type&&
               strcmp(surface_calpha_atom[j].res_id,t_m_n.t.res_id)==0)
             {cout<<"it is self"<<endl;}
            else
             {    
             new_n.res_num=surface_calpha_atom[j].res_num;
             new_n.res_id=new char[5];  
             new_n.res_id=surface_calpha_atom[j].res_id;
             new_n.chain_type=surface_calpha_atom[j].chain_type;
             new_n.x=surface_calpha_atom[j].x;
             new_n.y=surface_calpha_atom[j].y;
             new_n.z=surface_calpha_atom[j].z;

             new_n.d=sqrt(pow((surface_calpha_atom[i].x-surface_calpha_atom[j].x),2)+
                          pow((surface_calpha_atom[i].y-surface_calpha_atom[j].y),2)+
                          pow((surface_calpha_atom[i].z-surface_calpha_atom[j].z),2));
             t_m_n.n.push_back(new_n);
             }//
            }
            t_m_n_chain.push_back(t_m_n);
            t_m_n.n.clear(); 
         }

        for(i=0;i<t_m_n_chain.size();i++)
        {
         sort(t_m_n_chain[i].n, t_m_n_chain[i].n.size());
        }
        
         int p1,p2;
         double m,tx,ty,tz,ax,ay,az;


	for(p1=0;p1<t_m_n_chain.size();p1++)
           {
              //get geometrical center side chain of target residue
              for(p2=0;p2<w_c.size();p2++)
              {
               if(w_c[p2].res_num==t_m_n_chain[p1].t.res_num&&
                  w_c[p2].chain_type==t_m_n_chain[p1].t.chain_type&&
                  strcmp(w_c[p2].res_id,t_m_n_chain[p1].t.res_id)==0)
                {
                     tx=0;ty=0;tz=0;ax=0;ay=0;az=0;
                     m=0;
                     //aa_key=p2;

                  for(j=0;j<w_c[p2].atom_str.size();j++)
                  {
                     //cout<<w_c[p2].atom_str[j].atom_id<<endl;
                    if(strncmp(w_c[p2].atom_str[j].atom_id,"CA ",3)!=0&&strncmp(w_c[p2].atom_str[j].atom_id,"C  ",3)!=0
                     &&strncmp(w_c[p2].atom_str[j].atom_id,"N  ",3)!=0&&strncmp(w_c[p2].atom_str[j].atom_id,"O  ",3)!=0)
                     {
                       tx=tx+w_c[p2].atom_str[j].x;
                       ty=ty+w_c[p2].atom_str[j].y;
                       tz=tz+w_c[p2].atom_str[j].z;
                       m=m+1;
                      } //end if(strncmp
                  }// end for(j=0

                  if(m!=0)
                  {ax=tx/m;
                   ay=ty/m;
                   az=tz/m;
                   }

                } //end if(w_c[p2]
             }//end for(p2=0..

           t_m_n_chain[p1].t.sc_x=ax;
           t_m_n_chain[p1].t.sc_y=ay;
           t_m_n_chain[p1].t.sc_z=az;

        }//end for(p1=0



   

           Vector3D Ca;      //calpha atom 
           Vector3D GcSc;    //geometrical center of side chain
           Vector3D ScCa;     //Ca-GCSc

           Vector3D Pn_3;        // normal vector of triangle by FnAtom,SnAtom,Atom_3rd
        
           double x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D,global_angle;

	for(i=0;i<t_m_n_chain.size();i++)
           {
              Ca.SetXComponent(t_m_n_chain[i].t.x);
              Ca.SetYComponent(t_m_n_chain[i].t.y);
              Ca.SetZComponent(t_m_n_chain[i].t.z);
              GcSc.SetXComponent(t_m_n_chain[i].t.sc_x);
              GcSc.SetYComponent(t_m_n_chain[i].t.sc_y);
              GcSc.SetZComponent(t_m_n_chain[i].t.sc_z);
              ScCa=Ca-GcSc;

              x1=t_m_n_chain[i].n[0].x;
              y1=t_m_n_chain[i].n[0].y;
              z1=t_m_n_chain[i].n[0].z;

              x2=t_m_n_chain[i].n[1].x;
              y2=t_m_n_chain[i].n[1].y;
              z2=t_m_n_chain[i].n[1].z;

              x3=t_m_n_chain[i].n[2].x;
              y3=t_m_n_chain[i].n[2].y;
              z3=t_m_n_chain[i].n[2].z;

              A=y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2);
              B=z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2);
              C=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
              D=(-1)*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1)); 
  
              Pn_3.SetXComponent(A);
              Pn_3.SetYComponent(B);
              Pn_3.SetZComponent(C);
             
	     //angle between w_t_c and Plane_Normal
             global_angle=angle_between_two_vectors(ScCa,Pn_3);
             t_m_n_chain[i].t.angle=global_angle;    

          }



        for(i=0;i<t_m_n_chain.size();i++)
        {
          cout<<"target"<<" "<<t_m_n_chain[i].t.res_num<<" "<<t_m_n_chain[i].t.res_id<<" "
              <<t_m_n_chain[i].t.chain_type<<" "
              <<t_m_n_chain[i].t.x<<" "
              <<t_m_n_chain[i].t.y<<" "
              <<t_m_n_chain[i].t.z<<" "
              <<t_m_n_chain[i].t.sc_x<<" "
              <<t_m_n_chain[i].t.sc_y<<" "
              <<t_m_n_chain[i].t.sc_z<<" "
              <<t_m_n_chain[i].t.angle<<endl;

           /*
          for(j=0;j<t_m_n_chain[i].n.size();j++)
          {
           cout<<t_m_n_chain[i].n[j].res_num<<" "<<t_m_n_chain[i].n[j].res_id<<" "
               <<t_m_n_chain[i].n[j].chain_type<<" "
               <<t_m_n_chain[i].n[j].x<<" "
               <<t_m_n_chain[i].n[j].y<<" "
               <<t_m_n_chain[i].n[j].z<<endl;
          }
           */

        }


} 
