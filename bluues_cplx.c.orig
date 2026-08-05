#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

struct Area {
              double sas, ctc, pc_sas, pc_ctc, sas_self, ctc_self, sas_res, ctc_res;
            } Area;

struct Summary {
              double ccpc, pppc, ccsum, ppsum; 
              double sas, ses;
              double coul_nrg, solv_nrg, born_nrg, coul_solv_nrg, sasa_nrg, total_nrg;  
              double charge, charge_mol_2, charge_mol_1; 
            } Summary;

struct System_area {
               struct Area system;
               struct Area *atoms;
               struct Area *residues;
               struct Area *chains;
               struct Area *segments;
               struct Area *models;
             } System_area;

struct At_srf
{
int n_isrf; 
int *isrf;
};

struct Srf_pt
{
double *r;
double *vec;
double a;
double pot;
int at, atn;
int type;
};

struct Srf
{
int alloc_srf_pt;
int n_srf_pt_tot;
struct Srf_pt *srf_pt;
struct At_srf *at_srf;
double area;
};

struct System { 
               int n_atoms;
               struct Atom *atoms;
               int n_residues;
               struct Residue *residues;
               int n_chains;
               struct Chain *chains;
               int n_segments;
               struct Segment *segments;
               int n_models;
               struct Model *models;
              }  System;

struct Atom_grid {
                  int **atom_node;
                  int PBC_ON; 
                  double pbcx, pbcy, pbcz;
                  int *n_atom_node;
                  int grid_X, grid_Y, grid_Z, grid_size, max_atom_node;
                  double x_min, y_min, z_min, x_max, y_max, z_max;
                  double mesh;
                  double dx,dy,dz;
                 } Atom_grid;

struct Neighbour {
                  int max_n;             
                  int n_neighbours;             
                  int *list; 
                  double *d;
                 } Neighbour; 


struct Atom {
		char at_name[5]; 
		char alt_loc;
		char res_name[4];
		char chain; 
                char element[3]; 
		int model;
                int at_n;
		int res_n;
		char res_ins; 
                double coor[3]; 
                double ref[3];
		double occ,temp;
		char segid[5];
     		char pdb_chrg[3]; 
    		double mass;
    		char chem_group[5];  
    		char group[5];  
    		double radius;
    		double charge;
    		int atom_type;
		int i;
	    } Atom;

struct Trj {
           int nf, naxf;  
           double **coor; 
           } Trj;

struct Residue { char res_name[4];
                 char chain;
                 char res_ins;
                 int res_n;
                 int model;
                 char seg_name[5];
                 int beg, end ; 
                 int n_alt_loc;  
                 int prev, next;
                 double CA[3], CG[3], CM[3]; 
                 int res_type;
               } Residue;

struct Chain { char ch_name;
               int beg, end;
               int model;
               char seg_name[5];
             } Chain;


struct Segment { char seg_name[5];
               int beg, end;
               int model;
               }  Segment;

struct Model { int model_n;
               int beg, end;
               }  Model;

FILE *file_open(char *fname, char *acc);

void diffv(double *v, double *r2, double *r1);
void sumv(double *v, double *r2, double *r1);
void cmul(double *v, double a, double *v0);
void scale(double *v, double a);
double distv(double *r1, double *r2);
double distv_pbc(double *r1, double *r2, double *s);
double dist(double x1, double y1, double z1, double x2, double y2, double z2);

double anglev(double *x1, double *x2, double *x3);
double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

double torsionv(double *x1, double *x2, double *x3, double *x4);
double torsion(double x1, double y1, double z1, 
																				double x2, double y2, double z2,
                    double x3, double y3, double z3, 
	                   double x4, double y4, double z4);
	                    
double dotv(double *x1, double *x2);
double dot(double x1, double y1, double z1, double x2, double y2, double z2);
double dotvn(double *x1, double *x2,int n);
double modv(double *x1);
double modvn(double *x1, int n);

void vecv(double *x, double *x1, double *x2);
void vec(double *x, double *y, double *z,
         double x1, double y1, double z1, double x2, double y2, double z2);

int is_within(struct Atom *atoms,int n_atoms,struct Atom probe_atom, double cutoff);

void read_atom_pdb(char *buf, struct Atom *atom);

void read_PQR_atoms(FILE *fp1, int *n_atoms, struct Atom **atoms, struct Trj *trj);

void read_atom_pqr(char *buf, struct Atom *atom);

void write_PDB_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj);

void write_PQR_atoms(FILE *fp1, int n_atoms, struct Atom *atoms, struct Trj trj);

void write_atom_pdb(char *buf, struct Atom atom);

void write_atom_pqr(char *buf, struct Atom atom);

void copy_atom(struct Atom *atom_to, struct Atom atom_from);

void get_grid_parameters(struct System sistema, double probe_radius,
                double cutoff, struct Atom_grid *atom_grid, char type);

void print_grid_info(struct Atom_grid atom_grid);

void atoms_on_grid(struct Atom_grid *atom_grid, struct System sistema);

void grid_neighbour(struct Atom_grid atom_grid, 
                struct System sistema, double cutoff, struct Neighbour *neighbour, char type);

void system_area_alloc(struct System system, struct System_area *system_area);

void system_area_free(struct System_area *system_area);

struct Neighbour *neighbour_alloc(int n, int max_nei);

void neighbour_realloc(struct Neighbour *neighbours, int n);

void neighbour_free(struct Neighbour *neighbours, int n);

void atom_grid_alloc(struct Atom_grid *atom_grid);

void atom_grid_realloc(struct Atom_grid *atom_grid);

void atom_grid_free(struct Atom_grid *atom_grid);

void srf_free(struct Srf *srf, int n_atoms);

void cmul(double *v, double a, double *v0)
{
	v[0] = v0[0] * a;
	v[1] = v0[1] * a;
	v[2] = v0[2] * a;
}

void scale(double *v, double a)
{
	v[0] = v[0] * a;
	v[1] = v[1] * a;
	v[2] = v[2] * a;
}

void diffv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] - r1[0];
	v[1] = r2[1] - r1[1];
	v[2] = r2[2] - r1[2];
}

void sumv(double *v, double *r2, double *r1)
{
	v[0] = r2[0] + r1[0];
	v[1] = r2[1] + r1[1];
	v[2] = r2[2] + r1[2];
}

double dotv(double *x1, double *x2)
{       
    double d;
    d = x1[0]*x2[0] +  x1[1]*x2[1] + x1[2]*x2[2];
    return d;
}

double modv(double *x1)
{       
    double d;
    d = sqrt(x1[0]*x1[0] +  x1[1]*x1[1] + x1[2]*x1[2]);
    return d;
}

double dotvn(double *x1, double *x2, int n)
{       
    int i;
    double d=0.0;
    for(i=0; i<n; i++)
    d = d + x1[i]*x2[i];
    return d;
}

double modvn(double *x1, int n)
{       
    int i;
    double d=0.0;
    for(i=0; i<n; i++)
    d = d + x1[i]*x1[i];
    d = sqrt(d);
    return d;
}

void vecv(double *x, double *x1, double *x2)
{
    x[0] = x1[1]*x2[2]-x1[2]*x2[1];
    x[1] = x1[2]*x2[0]-x1[0]*x2[2];
    x[2] = x1[0]*x2[1]-x1[1]*x2[0];
}
 
double distv(double *r1, double *r2)
{
	double d,v[3];
	diffv(v,r2,r1);
	d = dotv(v, v );
	d =  sqrt(d);
	return d;
}

double distv_pbc(double *r1, double *r2, double *s)
{
        int i; 
	double d,v[3];
	diffv(v,r2,r1);
        for(i = 0; i<3; i++)
          if(v[i] > s[i] / 2) v[i] = v[i] - s[i];
          else
          if(v[i] < -s[i] / 2) v[i] = v[i] + s[i];

	d = dotv(v, v );
	d =  sqrt(d);
	return d;
}

double anglev(double *x1, double *x2, double *x3)
{
      
    double d21,d32,x21[3],x32[3], f, angle;
    diffv(x21, x2, x1);
    diffv(x32, x3, x2);

    d21 = sqrt(dotv(x21,x21));
    d32 = sqrt(dotv(x32,x32));
    f = -dotv(x21,x32)/(d21 * d32);

    angle = 180.0 * acos(f) / M_PI;
    return angle;
}

double torsionv(double *x1, double *x2, double *x3, double *x4)
{
        double d21,d32,d43,x21[3],x32[3],x43[3],bxc[3];
        double  ac,ab,bc,abxc, t;
        int i;


        diffv(x21, x2, x1);
        diffv(x32, x3, x2);
        diffv(x43, x4, x3);
        d21 = sqrt(dotv(x21,x21));
        d32 = sqrt(dotv(x32,x32));
        d43 = sqrt(dotv(x43,x43));

        for(i = 0; i < 3; i++)
        {
         x21[i] = x21[i]/d21;
         x32[i] = x32[i]/d32;
         x43[i] = x43[i]/d43;
        }

        ab = dotv(x21,x32);
        bc = dotv(x32,x43);
        ac = dotv(x21,x43);

        vecv(bxc, x32, x43);
        abxc=dotv(x21, bxc);

        t = 180.0 * atan2f(abxc, -ac + ab*bc )/M_PI;
        return t;

}

double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double d;
	d = dot(x1-x2, y1-y2, z1-z2, x1-x2, y1-y2, z1-z2);
 	d = sqrt(d);
	return d; 
}

double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
    double d1,d2,angle;

    d1 = dot(x1-x2, y1-y2, z1-z2,x1-x2, y1-y2, z1-z2);
    d2 = dot(x3-x2, y3-y2, z3-z2,x3-x2, y3-y2, z3-z2);
    angle = dot((x1-x2)/sqrt(d1), (y1-y2)/sqrt(d1), (z1-z2)/sqrt(d1), (x3-x2)/sqrt(d2), (y3-y2)/sqrt(d2), (z3-z2)/sqrt(d2));

    angle = 180.0 * acos(angle) / M_PI;
    return angle;
}

double torsion(double x1, double y1, double z1, 
                            double x2, double y2, double z2,
                            double x3, double y3, double z3, 
                            double x4, double y4, double z4)
{
    double d1,d2,d3,p1_x,p1_y,p1_z,p2_x,p2_y,p2_z,p3_x,p3_y,p3_z;
    double  ac,ab,bc,abxc,bxc_x,bxc_y,bxc_z;
    double torsion;

    d1 = dot(x2-x1,y2-y1,z2-z1, x2-x1,y2-y1,z2-z1);
    d1 =  sqrt(d1);
    d2 = dot(x3-x2,y3-y2,z3-z2, x3-x2,y3-y2,z3-z2);
    d2 =  sqrt(d2);
    d3 = dot(x4-x3,y4-y3,z4-z3, x4-x3,y4-y3,z4-z3);
    d3 =  sqrt(d3);

    p1_x = (x2-x1)/d1; p1_y = (y2-y1)/d1; p1_z = (z2-z1)/d1;
    p2_x = (x3-x2)/d2; p2_y = (y3-y2)/d2; p2_z = (z3-z2)/d2;
    p3_x = (x4-x3)/d3; p3_y = (y4-y3)/d3; p3_z = (z4-z3)/d3;

    ab = dot(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z);
    bc = dot(p2_x, p2_y, p2_z, p3_x, p3_y, p3_z);
    ac = dot(p1_x, p1_y, p1_z, p3_x, p3_y, p3_z);

    vec(&bxc_x, &bxc_y, &bxc_z, p2_x, p2_y, p2_z, p3_x, p3_y, p3_z);
    abxc = dot(p1_x, p1_y, p1_z, bxc_x, bxc_y, bxc_z);

    torsion = 180.0 * atan2f(abxc, -ac + ab*bc )/M_PI;
    return torsion;
}

double dot(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double d;
    d = x1*x2 +  y1*y2 + z1*z2;
    return d;
}

void vec(double *x, double *y, double *z,
         double x1, double y1, double z1, double x2, double y2, double z2)
{
    *x =  (y1 * z2 - z1 * y2);
    *y = -(x1 * z2 - z1 * x2);
    *z =  (x1 * y2 - y1 * x2);
}

void write_PDB_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj)
{
	int i,j,k,l;
	char buf[120];

        if(trj.naxf != n_atoms)
          {
          printf("the number of atoms per frame (%i) is not equal to number of atoms (%i)... \n", trj.naxf, n_atoms); 

          }
        for(j = 1,l=0; j<= trj.nf; j++)
        {
	fprintf(fp, "MODEL %8i\n", j);
	for(i=0;i< n_atoms; i++)
	{ 
                for(k = 0; k< 3; k++)
                atoms[i].coor[k] = trj.coor[l][k];
		write_atom_pdb(buf, atoms[i]);
		fprintf(fp, "%s", buf);
                l++;
	}
	fprintf(fp, "ENDMDL\n");
        }
}

void write_PQR_atoms(FILE *fp, int n_atoms, struct Atom *atoms, struct Trj trj)
{
	int i,j,k,l;
	char buf[120];

        if(trj.naxf != n_atoms)
          {
          printf("the number of atoms per frame (%i) is not equal to number of atoms (%i)... exiting\n", trj.naxf, n_atoms); 
          exit(0);
          }
        for(j = 1,l=0; j<= trj.nf; j++)
        {
	fprintf(fp, "MODEL %8i\n", j);
	for(i=0;i< n_atoms; i++)
	{ 
                for(k = 0; k< 3; k++)
                atoms[i].coor[k] = trj.coor[l][k];
		write_atom_pqr(buf, atoms[i]);
		fprintf(fp, "%s", buf);
                l++;
	}
	fprintf(fp, "ENDMDL\n");
        }
}


void write_atom_pdb(char *buf, struct Atom atom)
{
	int i;
	char temp[5] = "    ";
	if(strlen(atom.at_name) < 4)
		for(i = 1; i<= strlen(atom.at_name); i++) 
			temp[i] = atom.at_name[i-1]; 
	else strcpy(temp,atom.at_name);
	sprintf(buf,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2lf%6.2lf      %4s%2s%2s\n",
             atom.at_n%100000, temp, atom.alt_loc, atom.res_name, atom.chain,
             atom.res_n, atom.res_ins,atom.coor[0], atom.coor[1], atom.coor[2], 
             atom.occ, atom.temp, atom.segid, atom.element, atom.pdb_chrg); 
}

void write_atom_pqr(char *buf, struct Atom atom)
{
	int i;
	char temp[5] = "    ";
	if(strlen(atom.at_name) < 4)
		for(i = 1; i<= strlen(atom.at_name); i++) 
			temp[i] = atom.at_name[i-1]; 
	else strcpy(temp,atom.at_name);
	sprintf(buf,"ATOM  %5i %4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f %7.4lf %7.4lf      %4s%2s%2s\n",
             atom.at_n, temp, atom.alt_loc, atom.res_name, atom.chain,
             atom.res_n, atom.res_ins,atom.coor[0], atom.coor[1], atom.coor[2], 
             atom.charge, atom.radius, atom.segid, atom.element, atom.pdb_chrg); 
}

void copy_atom(struct Atom *atom_to, struct Atom atom_from)
{
int i;
strcpy((*atom_to).at_name,atom_from.at_name);
(*atom_to).alt_loc = atom_from.alt_loc;
strcpy((*atom_to).res_name,atom_from.res_name);
(*atom_to).chain = atom_from.chain;
strcpy((*atom_to).element,atom_from.element);
(*atom_to).model = atom_from.model;
(*atom_to).at_n = atom_from.at_n;
(*atom_to).res_n = atom_from.res_n;
(*atom_to).res_ins = atom_from.res_ins;

for(i=0;i<3;i++)
{

(*atom_to).coor[i] = atom_from.coor[i];

}
(*atom_to).occ = atom_from.occ;
(*atom_to).temp = atom_from.temp;
strcpy((*atom_to).segid,atom_from.segid);
strcpy((*atom_to).pdb_chrg,atom_from.pdb_chrg);
(*atom_to).mass = atom_from.mass;
strcpy((*atom_to).chem_group,atom_from.chem_group);
strcpy((*atom_to).group,atom_from.group);
(*atom_to).radius = atom_from.radius;
(*atom_to).charge = atom_from.charge;
(*atom_to).atom_type = atom_from.atom_type;
(*atom_to).i = atom_from.i;
}

void read_PQR_atoms(FILE *fp1, int *n_atoms, struct Atom *(*atoms), struct Trj *trj)
{
	char buf[120];
	int i=0, k, n_models,is_trj;
        int mod_id=0;
        struct Atom tmp_atom;


	while(fgets(buf,120,fp1) != NULL )
        {
    	if(!strncmp("ATOM",buf,4) && mod_id==0) i++;
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
        }
*n_atoms = i; 
n_models = mod_id;
if(mod_id <= 1) mod_id = 1;
  (*atoms) = (struct Atom *) calloc((size_t) *n_atoms , sizeof(struct Atom));
  if((*atoms) == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", *n_atoms);
     exit(0);
    }
  (*trj).nf = mod_id; 
  (*trj).naxf = *n_atoms; 
  trj->coor = (double **) calloc((size_t) ((*trj).nf * (*trj).naxf), sizeof(double *));
  if(trj->coor == NULL) 
    {
     printf("could not allocate memory for %i atoms... exiting...\n", (*trj).nf * (*trj).naxf);
     exit(0);
    }
   for(i = 0; i < (*trj).nf * (*trj).naxf; i++)
    {
    trj->coor[i] =  (double *) calloc((size_t) 3, sizeof(double));
    if(trj->coor[i] == NULL) 
    {
     printf("could not allocate memory for %i-th atom coordinates... exiting...\n", i);
     exit(0);
    }
    }
        rewind(fp1);
        i = 0;
	while(fgets(buf,120,fp1) != NULL)
    	if(!strncmp("ATOM",buf,4)) 
            {
                        if(i < *n_atoms)
                        {
			read_atom_pqr(buf, &((*atoms)[i]));
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = (*atoms)[i].coor[k];
                        i++;
                        }
                        else if(mod_id > 1)
                        {
			read_atom_pqr(buf, &tmp_atom);
                        for(k=0;k<3;k++)
                        (*trj).coor[i][k] = tmp_atom.coor[k];
                        i++;
                        }
            if(!(i%10000))  printf("%i atoms read\n",i);
	    }
	    else
	      if(!strncmp("ENDMDL",buf,6)) 
              mod_id++; 
}


void read_atom_pqr(char *buf, struct Atom *atom)
{

	char at_rec[5];
    char tok[10];              
 
    strncpy(tok,buf,4);
    tok[4] = '\0';
    sscanf(tok,"%s", at_rec);
    if(strncmp("ATOM",at_rec,4))
    {
    	printf("The ATOM line does not start with string ATOM... exiting...\n");
        exit(1);
    }

    strncpy(tok,buf + 6,5);
    tok[5] = '\0';
    sscanf(tok,"%i",&(atom->at_n));

    strncpy(tok,buf + 12,4);
    tok[4] = '\0';
    sscanf(tok,"%s", atom->at_name);
 
    strncpy(tok,buf + 16,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->alt_loc)) == -1) atom->alt_loc=' ';
 

	strncpy(tok,buf + 17,3);
    tok[3] = '\0'; 
    sscanf(tok,"%s", atom->res_name);

	strncpy(tok,buf + 21,1);
    tok[1] = '\0';
    if(sscanf(tok,"%c", &(atom->chain)) == EOF) atom->chain = ' ';

    strncpy(tok,buf + 22,4);
    tok[4] = '\0';
    sscanf(tok,"%i", &(atom->res_n));

	strncpy(tok,buf + 26,1);
    tok[1] = '\0';
    if (sscanf(tok,"%c", &(atom->res_ins)) == EOF) atom->res_ins=' ';

    strncpy(tok,buf + 30,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[0]));

	strncpy(tok,buf + 38,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[1]));

    strncpy(tok,buf + 46,8);
    tok[8] = '\0';
    sscanf(tok,"%lf", &(atom->coor[2]));

    sscanf(buf + 55 ,"%lf %lf", &(atom->charge), &(atom->radius));





}


struct Neighbour *neighbour_alloc(int n, int max_nei)
{
    int i;
    struct Neighbour *p;
    p = (struct Neighbour *) calloc((size_t) n , sizeof(struct Neighbour));
    if(p==NULL) printf("Could not allocate memory for Neighbour... Exiting...\n");

    for(i = 0; i< n; i++)
    p[i].max_n = max_nei;

    for(i = 0; i< n; i++)
    {
    p[i].list = (int *) calloc((size_t) max_nei, sizeof(int));
    if(p[i].list == NULL) 
    { printf("Could not allocate memory Neighbour.list[%i]... Exiting...\n", i);    exit(0);}
    }
    for(i = 0; i< n; i++)
    {
    p[i].d = (double *) calloc((size_t) max_nei, sizeof(double));
    if(p[i].d == NULL) 
    { printf("Could not allocate memory Neighbour.d[%i]... Exiting...\n", i);    exit(0);}
    }
    return p;
}

void neighbour_free(struct Neighbour *neighbours, int n)
{
    int i;
    for(i = 0; i< n; i++)
    free((neighbours)[i].list);
    free(neighbours);
}

void neighbour_realloc(struct Neighbour *neighbours, int n)
{
    int i;
    for(i = 0; i< n; i++)
    {
    neighbours[i].max_n = neighbours[i].n_neighbours;
    neighbours[i].list = (int *) realloc(neighbours[i].list, (size_t) neighbours[i].n_neighbours * sizeof(int));
    neighbours[i].d = (double *) realloc(neighbours[i].d, (size_t) neighbours[i].n_neighbours * sizeof(double));
    }
}


void atom_grid_free(struct Atom_grid *atom_grid)
{
int i;
for(i=0; i< (*atom_grid).grid_size; i++)
    free((*atom_grid).atom_node[i]);
free((*atom_grid).atom_node); 
free((*atom_grid).n_atom_node); 
}

void atom_grid_alloc(struct Atom_grid *atom_grid)
{
int i;

atom_grid->atom_node = (int **) calloc((size_t) (*atom_grid).grid_size , sizeof(int*));
if(atom_grid->atom_node == NULL)
    { printf("Could not allocate memory for atom_grid->atom_node... Exiting...\n");    exit(0);}

atom_grid->n_atom_node = (int *) calloc((size_t) (*atom_grid).grid_size , sizeof (int ));
if(atom_grid->n_atom_node == NULL)
    { printf("Could not allocate memory for atom_grid->n_atom_node... Exiting...\n");    exit(0);}
for(i = 0; i< (*atom_grid).grid_size; i++)
 {

(*atom_grid).atom_node[i] = (int *) calloc((size_t) (*atom_grid).max_atom_node , sizeof(int));
if(atom_grid->atom_node[i] == NULL)
    { printf("Could not allocate memory for atom_grid->atom_node[%i]... Exiting...\n", i);    exit(0);}
 }
atom_grid->PBC_ON = 0;
}

void atom_grid_realloc(struct Atom_grid *atom_grid)
{
int i;
for(i=0; i< (*atom_grid).grid_size;i++)
atom_grid->atom_node[i] = (int *) realloc(atom_grid->atom_node[i], (size_t) (*atom_grid).n_atom_node[i] * sizeof(int));
}

FILE *file_open(char *fname,char *acc) {
    FILE *fp;
    fp =fopen(fname,acc);
    if (fp==NULL) 
	{
        fprintf(stderr,"unable to open file %s\n",fname);
        exit(1);
    }
    return(fp);
}

struct Flag_par {  
	       char file_in_pdb[256];
	       char file_out[256];
	       char file_msms[256];
	       char file_ns[256];
               char chain_p[256], chain_l[256];
               double msms_area;
               double probe_radius;
               double test_charge_radius;
	       double inflate;
               double min_radius;
               double salt_radius;
               double pdie;
               double sdie;
               double ions;
               double temp;
               double cutoff;
               double cutoff2;
               double cutoff_short;
               double maxf;
               double gamma;
               double kp;
               int gbr;
               int surface;
               int srfpot;
               int print_surface;
               int print_srfpot;
               int srfpotav;
               int solv_e;
               int pqg;
               int grid, nx, ny, nz;
	       double mesh;
               int msms;
               int ns;
               int rpns;
               int verbose;
	     } Flag_par;

struct Area_atom {
        double ses, sas;
        double area;
        int n_model;
        char  segment[5]; 
        char chain; 
	int res_n;
	char res_ins; 
        char  at_name[5]; 
        char alt_loc;
        char res_name[4];
        int at_n;
        } Area_atom;

void grid_srf_dens(
                struct System system,
                struct Atom_grid atom_grid,
                struct Flag_par flag_par,
                struct Srf *srf);

struct Msms_face
{
int n1,n2,n3,type;
};

struct Energy
{
double area, solv, coul, pol, born, tot; 
double *area_a, *coul_a, *pol_a, *solv_a, *born_a, *tot_a;
};

struct Const
{
double q, N_av, e0, kb, RT, kd, k_el;
};

void get_srf(struct System sys, struct Atom_grid *patom_grid, struct Srf *psrf_ses, struct Srf *psrf_sas, struct Flag_par flag_par);

	
void pqr2gbr6(struct System sys, struct Srf srf_ses, struct Srf Srf_sas, struct Flag_par flag_par, struct Neighbour *neighbours, double *gbr6_ses, double *gbr6_sas);

void print_radii(struct System sys, struct Trj trj, double *gbr6_ses, double *gbr6_sas, struct Flag_par flag_par);

void print_srf_pot(struct System sys, struct Trj trj, double *gbr6, struct Srf *srf_ses, struct Srf *srf_sas, struct Neighbour *neighbours, struct Area_atom *area_atom, struct Flag_par flag_par);

void print_summary(struct Summary sum, struct Flag_par flag_par);


void srfpot(struct System *sys, struct Srf *srf_ses, struct Srf *srf, double *gbr6, struct Neighbour *neighbours, struct Flag_par flag_par);

void calc_nrg(struct System sys, struct Energy *pnrg, double *gbr6, struct Srf srf, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par);

void print_nrg(struct Energy nrg, struct System sys, struct Flag_par flag_par);

void out_grid(struct System sys, struct Atom_grid atom_grid, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par);

void scalev(double t, double *v);

double fgb(double rgb1, double rgb2, double rw, double rs, double d, 
           double pdie, double sdie, double kd, double kp);

int cmp_area_atom_by_string(const void *p1, const void *p2);
int cmp_area_atom_by_n(const void *p1, const void *p2);
int cmp_int(const void *p1, const void *p2);

void init_flag_par(struct Flag_par *flag_par);
void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par);
void print_info_flag_par(struct Flag_par flag_par);




int main(int argc, char *argv[]) {
      FILE *fp, *fp1, *fp2, *fp3, *fp4;
      struct Flag_par flag_par;
      struct Const constants;
      struct Summary summary;
      struct System sys, mol_1, mol_2;
      struct Trj trj, trj_mol_1, trj_mol_2;
      struct Atom_grid atom_grid, atom_grid_mol_1, atom_grid_mol_2;
      struct Neighbour *neighbours_mol_2, *neighbours_mol_1,  *neighbours; 
      struct Energy nrg, nrg_mol_2, nrg_mol_1, nrg_diff;
      int i,j,k,l,m,n,h,i1,j1,i2,j2, *is_in_int_1, *is_in_int_2; 
      double d;
      double *gbr6_ses, *gbr6_sas;
      double *gbr6_ses_mol_2, *gbr6_sas_mol_2;
      double *gbr6_ses_mol_1, *gbr6_sas_mol_1;
      struct Srf *psrf_ses, *psrf_sas;
      struct Srf srf_ses, srf_ses_mol_2, srf_ses_mol_1;
      struct Srf srf_sas, srf_sas_mol_2, srf_sas_mol_1;
      double charge, charge_mol_2, charge_mol_1; 
      struct Area_atom *area_atom, *area_atom_mol_2, *area_atom_mol_1;
      char   buf[1048], buf_ns[4096];
      double Sxy, Sx, Sy, Sx2, Sy2, nn;
      double Sxyc, Sxc, Syc, Sx2c, Sy2c, nnc;
char *ch="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

        
        init_flag_par(&flag_par);      
        
        check_cmd_line(argc, argv, &flag_par);      
        
        print_info_flag_par(flag_par);      

      constants.q = 1.602176487E-19;
      constants.N_av = 6.02214179E23;
      constants.e0 = 8.8541878176E-12;
      constants.kb = 1.3806504E-23;
      constants.RT = 1.3806504 * 6.02214179 * flag_par.temp;
      constants.kd = 1E-10 * sqrt((1000 * 2 * flag_par.ions * constants.N_av * constants.q * constants.q )/ 
                         (flag_par.sdie * constants.e0 * constants.kb * flag_par.temp));
      constants.k_el = constants.N_av * 1E10 * constants.q * constants.q /(4 * M_PI * constants.e0);
      printf("k_el = %e (kJ/mol) %e (kcal/mol)\n", 
      constants.k_el/1000, constants.k_el/4184); 
      printf("Debye length = %e (A)\n", 1.0/constants.kd);

        fp = file_open(flag_par.file_in_pdb,"r");
        read_PQR_atoms(fp, &(sys.n_atoms), &(sys.atoms), &trj);
        fclose(fp);
        printf("\n%i total atoms read\n", sys.n_atoms);

        for(i = 0; i < sys.n_atoms; i++)
        for(k = 0; k < 3; k++)
        sys.atoms[i].coor[k]  = trj.coor[i][k];
        
        for(i=0; i< sys.n_atoms; i++)
        if(sys.atoms[i].radius < flag_par.min_radius) sys.atoms[i].radius = flag_par.min_radius;
        
        for(i=0; i< sys.n_atoms; i++)
        sys.atoms[i].radius = sys.atoms[i].radius + flag_par.inflate;

//        if(flag_par.verbose)
//        for(i=0; i < sys.n_atoms; i++) 
//        printf("atom %i: charge %lf  radius %lf\n", i, sys.atoms[i].charge, sys.atoms[i].radius);


        for(i=0,k=0; i< sys.n_atoms; i++)
        for(m=0; m < strlen(flag_par.chain_l); m++)
        if(sys.atoms[i].chain == flag_par.chain_l[m]) k++;
        mol_2.n_atoms = k;
        mol_2.atoms = calloc(mol_2.n_atoms, sizeof( struct Atom)); 
        trj_mol_2.coor = calloc(mol_2.n_atoms, sizeof(double *));
        for(k = 0; k<mol_2.n_atoms; k++)
          trj_mol_2.coor[k] = calloc(3, sizeof(double));
        for(i=0,k=0; i< sys.n_atoms; i++)
        for(m=0; m < strlen(flag_par.chain_l); m++)
        if(sys.atoms[i].chain == flag_par.chain_l[m]) 
          {
        copy_atom(&(mol_2.atoms[k]),sys.atoms[i]);
        for(l = 0; l < 3; l++)
        mol_2.atoms[k].coor[l] = trj.coor[i][l];
        for(l = 0; l < 3; l++)
        trj_mol_2.coor[k][l]  = trj.coor[i][l];
        k++;
          }


        for(i=0,k=0; i< sys.n_atoms; i++)
        for(m=0; m < strlen(flag_par.chain_p); m++)
        if(sys.atoms[i].chain == flag_par.chain_p[m]) k++;
        mol_1.n_atoms = k;
        mol_1.atoms = calloc(mol_1.n_atoms, sizeof( struct Atom)); 
        trj_mol_1.coor = calloc(mol_1.n_atoms, sizeof(double *));
        for(k = 0; k<mol_1.n_atoms; k++)
          trj_mol_1.coor[k] = calloc(3, sizeof(double));

        for(i=0,k=0; i< sys.n_atoms; i++)
        for(m=0; m < strlen(flag_par.chain_p); m++)
        if(sys.atoms[i].chain == flag_par.chain_p[m])
        if(is_within(mol_2.atoms,mol_2.n_atoms,sys.atoms[i],2*flag_par.cutoff))
        {
        copy_atom(&(mol_1.atoms[k]),sys.atoms[i]);
        for(l = 0; l < 3; l++)
        mol_1.atoms[k].coor[l] = trj.coor[i][l];
        for(l = 0; l < 3; l++)
        trj_mol_1.coor[k][l]  = trj.coor[i][l];
        k++;
        }
        for(j = k; j<mol_1.n_atoms; j++)
          free(trj_mol_1.coor[j]);
        mol_1.n_atoms = k;
        mol_1.atoms = realloc(mol_1.atoms, mol_1.n_atoms * sizeof( struct Atom)); 
        trj_mol_1.coor = realloc(trj_mol_1.coor, mol_1.n_atoms * sizeof(double *));
        trj_mol_1.nf = 1;
        trj_mol_1.naxf = mol_1.n_atoms;

        for(i=0,k=0; i< sys.n_atoms; i++)
        for(m=0; m < strlen(flag_par.chain_l); m++)
        if(sys.atoms[i].chain == flag_par.chain_l[m])
        if(is_within(mol_1.atoms,mol_1.n_atoms,sys.atoms[i],2*flag_par.cutoff))
        {
        copy_atom(&(mol_2.atoms[k]),sys.atoms[i]);
        for(l = 0; l < 3; l++)
        mol_2.atoms[k].coor[l] = trj.coor[i][l];
        for(l = 0; l < 3; l++)
        trj_mol_2.coor[k][l]  = trj.coor[i][l];
        k++;
        }
        for(j = k; j<mol_2.n_atoms; j++)
          free(trj_mol_2.coor[j]);
        mol_2.n_atoms = k;
        mol_2.atoms = realloc(mol_2.atoms, mol_2.n_atoms * sizeof( struct Atom)); 
        trj_mol_2.coor = realloc(trj_mol_2.coor, mol_2.n_atoms * sizeof(double *));
        trj_mol_2.nf = 1;
        trj_mol_2.naxf = mol_2.n_atoms;


        for(k=0,i=0; k< mol_1.n_atoms; k++)
        {
        copy_atom(&(sys.atoms[i]),mol_1.atoms[k]);
        for(l=0;l<3; l++)
        trj.coor[i][l]  = trj_mol_1.coor[k][l];
        for(l=0;l<3; l++)
        sys.atoms[i].coor[l]  = trj_mol_1.coor[k][l];
        i++; 
        }

        for(k=0; k< mol_2.n_atoms; k++)
        {
        copy_atom(&(sys.atoms[i]),mol_2.atoms[k]);
        for(l=0;l<3; l++)
        trj.coor[i][l]  = trj_mol_2.coor[k][l];
        for(l=0;l<3; l++)
        sys.atoms[i].coor[l]  = trj_mol_2.coor[k][l];
        i++; 
        }
        for(j = i; j<sys.n_atoms; j++)
          free(trj.coor[j]);
        sys.n_atoms = i;

        sys.atoms = realloc(sys.atoms, sys.n_atoms * sizeof( struct Atom)); 
        trj.coor = realloc(trj.coor, sys.n_atoms * sizeof(double *));
        trj.naxf = sys.n_atoms;

atom_grid.PBC_ON = 0;
atom_grid_mol_1.PBC_ON = 0;
atom_grid_mol_2.PBC_ON = 0;
get_grid_parameters(sys, flag_par.probe_radius, flag_par.cutoff, &atom_grid, 'v');
get_grid_parameters(sys, flag_par.probe_radius, flag_par.cutoff, &atom_grid_mol_1, 'v');
get_grid_parameters(sys, flag_par.probe_radius, flag_par.cutoff, &atom_grid_mol_2, 'v');

atom_grid_alloc(&atom_grid);
atom_grid_alloc(&atom_grid_mol_1);
atom_grid_alloc(&atom_grid_mol_2);
printf("Initial grid allocation for %i atoms: %i nodes\n", sys.n_atoms, atom_grid.max_atom_node);

neighbours = neighbour_alloc(sys.n_atoms, atom_grid.max_atom_node);
neighbours_mol_2 = neighbour_alloc(mol_2.n_atoms, atom_grid.max_atom_node);
neighbours_mol_1 = neighbour_alloc(mol_1.n_atoms, atom_grid.max_atom_node);
        if(flag_par.srfpot)
        {
        area_atom=(struct Area_atom *) calloc((size_t) sys.n_atoms, sizeof(struct Area_atom));
        area_atom_mol_2=(struct Area_atom *) calloc((size_t) mol_2.n_atoms, sizeof(struct Area_atom));
        area_atom_mol_1=(struct Area_atom *) calloc((size_t) mol_1.n_atoms, sizeof(struct Area_atom));
        }
        if(flag_par.solv_e || 1)
        {
        nrg.solv_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.solv_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.solv_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.solv_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        nrg.area_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.area_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.area_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.area_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        nrg.coul_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.coul_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.coul_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.coul_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        nrg.pol_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.pol_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.pol_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.pol_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));

        nrg.born_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.born_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.born_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.born_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        nrg.tot_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_diff.tot_a=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        nrg_mol_2.tot_a=(double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        nrg_mol_1.tot_a=(double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        }
        
        gbr6_ses = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
        gbr6_sas = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
        gbr6_ses_mol_2 = (double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        gbr6_sas_mol_2 = (double *) calloc((size_t) mol_2.n_atoms, sizeof(double));
        gbr6_ses_mol_1 = (double *) calloc((size_t) mol_1.n_atoms, sizeof(double));
        gbr6_sas_mol_1 = (double *) calloc((size_t) mol_1.n_atoms, sizeof(double));


	
          charge=0.0;
          charge_mol_2=0.0;
          charge_mol_1=0.0;

          for(i=0; i< sys.n_atoms; i++)
               charge = charge + sys.atoms[i].charge;

          for(i=0; i< mol_1.n_atoms; i++)
               charge_mol_1 = charge_mol_1 + mol_1.atoms[i].charge;
          for(i=0; i< mol_2.n_atoms; i++)
               charge_mol_2 = charge_mol_2 + mol_2.atoms[i].charge;
summary.charge = charge;
summary.charge_mol_1 = charge_mol_1;
summary.charge_mol_2 = charge_mol_2;
{
if(flag_par.verbose)
{
strcpy(buf, flag_par.file_out);
strcat(buf,"_int.pqr");
fp=file_open(buf,"w");
write_PQR_atoms(fp, sys.n_atoms, sys.atoms, trj);
fclose(fp);
strcpy(buf, flag_par.file_out);
strcat(buf,"_int_mol_1.pqr");
fp=file_open(buf,"w");
write_PQR_atoms(fp, mol_1.n_atoms, mol_1.atoms, trj_mol_1);
fclose(fp);
strcpy(buf, flag_par.file_out);
strcat(buf,"_int_mol_2.pqr");
fp=file_open(buf,"w");
write_PQR_atoms(fp, mol_2.n_atoms, mol_2.atoms, trj_mol_2);
fclose(fp);
}
}


atoms_on_grid(&atom_grid, sys);
atom_grid_realloc(&atom_grid); 
grid_neighbour(atom_grid, sys, flag_par.cutoff, neighbours,'v');
neighbour_realloc(neighbours, sys.n_atoms);
psrf_ses= &srf_ses;
psrf_sas= &srf_sas;
get_srf(sys, &atom_grid, psrf_ses, psrf_sas, flag_par);
if(flag_par.ns)
{
sprintf(buf_ns,"cp %s.area %s_cplx.area; cp %s.face %s_cplx.face; cp %s.prm  %s_cplx.prm; cp %s.vert %s_cplx.vert; cp %s.xyzr %s_cplx.xyzr",
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns
);
if(system(buf_ns))
{ printf("the call: %s failed.... check it in your terminal\nexiting.....\n",buf_ns); exit(0);}
}
else if(flag_par.msms)
{
sprintf(buf_ns,"cp %s.face %s_cplx.face; cp %s.vert %s_cplx.vert; cp %s.xyzr %s_cplx.xyzr",
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms
);
if(system(buf_ns))
{ printf("the call: %s failed.... check it in your terminal\nexiting.....\n",buf_ns); exit(0);}
}

pqr2gbr6(sys, srf_ses, srf_sas, flag_par, neighbours, gbr6_ses, gbr6_sas);
srfpot(&sys, &srf_ses, &srf_sas, gbr6_ses, neighbours, flag_par); 
print_srf_pot(sys, trj, gbr6_ses, &srf_ses, &srf_sas, neighbours, area_atom, flag_par);
calc_nrg(sys, &nrg, gbr6_ses, srf_sas, neighbours, constants, flag_par);
if(flag_par.verbose)
print_radii(sys, trj, gbr6_ses, gbr6_sas, flag_par);
if(flag_par.verbose)
print_nrg(nrg,sys,flag_par);

strcpy(buf, flag_par.file_out);
strcat(flag_par.file_out,"_mol_1");
atoms_on_grid(&atom_grid_mol_1, mol_1);
atom_grid_realloc(&atom_grid_mol_1); 
grid_neighbour(atom_grid_mol_1, mol_1, flag_par.cutoff, neighbours_mol_1,'v');
neighbour_realloc(neighbours_mol_1, mol_1.n_atoms);
psrf_ses= &srf_ses_mol_1;
psrf_sas= &srf_sas_mol_1;
get_srf(mol_1, &atom_grid_mol_1, psrf_ses, psrf_sas, flag_par);
if(flag_par.ns)
{
sprintf(buf_ns,
"cp %s.area %s_mol_1.area; cp %s.face %s_mol_1.face; cp %s.prm  %s_mol_1.prm; cp %s.vert %s_mol_1.vert; cp %s.xyzr %s_mol_1.xyzr",
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns
);
if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}
else if(flag_par.msms)
{
sprintf(buf_ns,
"cp %s.face %s_mol_1.face; cp %s.vert %s_mol_1.vert; cp %s.xyzr %s_mol_1.xyzr",
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms
);
if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}

pqr2gbr6(mol_1, srf_ses_mol_1, srf_sas_mol_1, flag_par, neighbours_mol_1, gbr6_ses_mol_1, gbr6_sas_mol_1);
srfpot(&mol_1, &srf_ses_mol_1, &srf_sas_mol_1, gbr6_ses_mol_1, neighbours_mol_1, flag_par); 
print_srf_pot(mol_1, trj_mol_1, gbr6_ses_mol_1, &srf_ses_mol_1, &srf_sas_mol_1, neighbours_mol_1, area_atom_mol_1, flag_par);
calc_nrg(mol_1,&nrg_mol_1, gbr6_ses_mol_1, srf_sas_mol_1, neighbours_mol_1, constants, flag_par);
if(flag_par.verbose)
print_nrg(nrg_mol_1,mol_1,flag_par);
if(flag_par.verbose)
print_radii(mol_1, trj_mol_1, gbr6_ses_mol_1, gbr6_sas_mol_1, flag_par);

strcpy(flag_par.file_out,buf);
strcpy(buf, flag_par.file_out);
strcat(flag_par.file_out,"_mol_2");
atoms_on_grid(&atom_grid_mol_2, mol_2);
atom_grid_realloc(&atom_grid_mol_2); 
grid_neighbour(atom_grid_mol_2, mol_2, flag_par.cutoff, neighbours_mol_2,'v');
neighbour_realloc(neighbours_mol_2, mol_2.n_atoms);
psrf_ses= &srf_ses_mol_2;
psrf_sas= &srf_sas_mol_2;
get_srf(mol_2, &atom_grid_mol_2, psrf_ses, psrf_sas, flag_par);
if(flag_par.ns)
{
sprintf(buf_ns,"cp %s.area %s_mol_2.area; cp %s.face %s_mol_2.face; cp %s.prm  %s_mol_2.prm; cp %s.vert %s_mol_2.vert; cp %s.xyzr %s_mol_2.xyzr",
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns
);
if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}
else if(flag_par.msms)
{
sprintf(buf_ns,
"cp %s.face %s_mol_2.face; cp %s.vert %s_mol_2.vert; cp %s.xyzr %s_mol_2.xyzr",
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms
);
if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}

if(flag_par.ns)
{
sprintf(buf_ns,"mv %s_cplx.area %s.area; mv %s_cplx.face %s.face; mv %s_cplx.prm  %s.prm; mv %s_cplx.vert %s.vert; mv %s_cplx.xyzr %s.xyzr",
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns
);

if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}
else if(flag_par.msms)
{
sprintf(buf_ns,"mv %s_cplx.face %s.face; mv %s_cplx.vert %s.vert; mv %s_cplx.xyzr %s.xyzr",
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms
);
if(system(buf_ns))
{ printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n", buf_ns); exit(0);}
}

pqr2gbr6(mol_2, srf_ses_mol_2, srf_sas_mol_2, flag_par, neighbours_mol_2, gbr6_ses_mol_2, gbr6_sas_mol_2);
srfpot(&mol_2, &srf_ses_mol_2, &srf_sas_mol_2, gbr6_ses_mol_2, neighbours_mol_2, flag_par); 
print_srf_pot(mol_2, trj_mol_2, gbr6_ses_mol_2, &srf_ses_mol_2, &srf_sas_mol_2, neighbours_mol_2, area_atom_mol_2, flag_par);
calc_nrg(mol_2,&nrg_mol_2, gbr6_ses_mol_2, srf_sas_mol_2, neighbours_mol_2, constants, flag_par);
if(flag_par.verbose)
print_nrg(nrg_mol_2,mol_2,flag_par);
if(flag_par.verbose)
print_radii(mol_2, trj_mol_2, gbr6_ses_mol_2, gbr6_sas_mol_2, flag_par);
strcpy(flag_par.file_out,buf);
summary.sasa_nrg = nrg_diff.area = nrg.area - nrg_mol_1.area - nrg_mol_2.area;
summary.coul_nrg = nrg_diff.coul = nrg.coul - nrg_mol_1.coul - nrg_mol_2.coul;
summary.solv_nrg = nrg_diff.solv = nrg.solv - nrg_mol_1.solv - nrg_mol_2.solv;
summary.coul_solv_nrg = nrg_diff.pol = nrg.pol - nrg_mol_1.pol - nrg_mol_2.pol;
summary.born_nrg = nrg_diff.born = nrg.born - nrg_mol_1.born - nrg_mol_2.born;
summary.total_nrg = nrg_diff.tot = nrg.tot - nrg_mol_1.tot - nrg_mol_2.tot;
summary.sas = srf_sas.area - srf_sas_mol_1.area - srf_sas_mol_2.area;
summary.ses = srf_ses.area - srf_ses_mol_1.area - srf_ses_mol_2.area;
for(k=0,i=0; k<mol_1.n_atoms; k++)
{
nrg_diff.coul_a[i] = nrg.coul_a[i] - nrg_mol_1.coul_a[k];
nrg_diff.pol_a[i] = nrg.pol_a[i] - nrg_mol_1.pol_a[k];
nrg_diff.born_a[i] = nrg.born_a[i] - nrg_mol_1.born_a[k];
nrg_diff.area_a[i] = nrg.area_a[i] - nrg_mol_1.area_a[k];
nrg_diff.solv_a[i] = nrg.solv_a[i] - nrg_mol_1.solv_a[k];
nrg_diff.tot_a[i] = nrg.tot_a[i] - nrg_mol_1.tot_a[k];
i++;
}
for(k=0; k<mol_2.n_atoms; k++)
{
nrg_diff.coul_a[i] = nrg.coul_a[i] - nrg_mol_1.coul_a[k];
nrg_diff.pol_a[i] = nrg.pol_a[i] - nrg_mol_1.pol_a[k];
nrg_diff.born_a[i] = nrg.born_a[i] - nrg_mol_1.born_a[k];
nrg_diff.area_a[i] = nrg.area_a[i] - nrg_mol_1.area_a[k];
nrg_diff.solv_a[i] = nrg.solv_a[i] - nrg_mol_1.solv_a[k];
nrg_diff.tot_a[i] = nrg.tot_a[i] - nrg_mol_1.tot_a[k];
i++;
}
if(flag_par.verbose)
{
strcpy(buf, flag_par.file_out);
strcat(flag_par.file_out,"_diff");
print_nrg(nrg_diff,sys,flag_par);
strcpy(flag_par.file_out,buf);
}
grid_neighbour(atom_grid, sys, flag_par.cutoff_short, neighbours,'v');
neighbour_realloc(neighbours, sys.n_atoms);
strcpy(buf, flag_par.file_out);
strcat(buf, ".crg_pair");
if(flag_par.verbose)
fp1=file_open(buf,"w");
if(flag_par.verbose && (flag_par.msms || flag_par.ns))
{
strcpy(buf, flag_par.file_out);
strcat(buf, ".srf_pot_pair");
fp2=file_open(buf,"w");
strcpy(buf, flag_par.file_out);
strcat(buf, "_mol_1.srf");
fp3=file_open(buf,"w");
strcpy(buf, flag_par.file_out);
strcat(buf, "_mol_2.srf");
fp4=file_open(buf,"w");
}

is_in_int_1 = calloc(srf_ses_mol_1.n_srf_pt_tot, sizeof(int));
is_in_int_2 = calloc(srf_ses_mol_2.n_srf_pt_tot, sizeof(int));
nn=Sx=Sy=Sxy=Sx2=Sy2=0.0;
nnc=Sxc=Syc=Sxyc=Sx2c=Sy2c=0.0;

for(i=0,i1=0,i2=0,j1=-1,j2=-1; i< mol_1.n_atoms; i++)
{
for(j=0; j< neighbours[i].n_neighbours; j++) 
{
k = neighbours[i].list[j];
if(k >= mol_1.n_atoms)
{

if(flag_par.verbose)
fprintf(fp1, "%5i %4s %4s %4i %c crg: %7.3lf -- %5i %4s %4s %4i %c crg:  %7.3lf -- d: %10.3lf\n",
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               sys.atoms[i].charge,
               sys.atoms[k].at_n,
               sys.atoms[k].at_name,
               sys.atoms[k].res_name,
               sys.atoms[k].res_n,
               sys.atoms[k].chain,
               sys.atoms[k].charge,
               distv(sys.atoms[i].coor,sys.atoms[k].coor));
Sxc=Sxc + sys.atoms[i].charge;
Syc=Syc + sys.atoms[k].charge;
Sxyc=Sxyc + sys.atoms[i].charge * sys.atoms[k].charge;
Sx2c=Sx2c + sys.atoms[i].charge * sys.atoms[i].charge;
Sy2c=Sy2c + sys.atoms[k].charge * sys.atoms[k].charge;
nnc = nnc+1.0;

if(flag_par.msms || flag_par.ns)
for(l=0; l < srf_ses_mol_1.at_srf[i].n_isrf ; l++)
{
m =srf_ses_mol_1.at_srf[i].isrf[l];
for(n=0; n < srf_ses_mol_2.at_srf[k - mol_1.n_atoms].n_isrf ; n++)
{
h = srf_ses_mol_2.at_srf[k- mol_1.n_atoms].isrf[n];
d = distv( srf_ses_mol_1.srf_pt[m].r, srf_ses_mol_2.srf_pt[h].r);
if(d<= flag_par.cutoff_short)
{
if(flag_par.verbose)
fprintf(fp2, "%5i %4s %4s %4i %c pot: %7.3lf -- %5i %4s %4s %4i %c pot:  %7.3lf -- d: %10.3lf\n",
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               srf_ses_mol_1.srf_pt[m].pot,
               sys.atoms[k].at_n,
               sys.atoms[k].at_name,
               sys.atoms[k].res_name,
               sys.atoms[k].res_n,
               sys.atoms[k].chain,
               srf_ses_mol_2.srf_pt[h].pot,
               d);
is_in_int_1[m] = 1;
is_in_int_2[h] = 1;
Sx=Sx + srf_ses_mol_1.srf_pt[m].pot;
Sy=Sy + srf_ses_mol_2.srf_pt[h].pot;
Sxy=Sxy + srf_ses_mol_1.srf_pt[m].pot * srf_ses_mol_2.srf_pt[h].pot;
Sx2=Sx2 + srf_ses_mol_1.srf_pt[m].pot * srf_ses_mol_1.srf_pt[m].pot;
Sy2=Sy2 + srf_ses_mol_2.srf_pt[h].pot * srf_ses_mol_2.srf_pt[h].pot;
nn = nn+1.0;
}
} 
}
}
}
}

summary.ccpc = (nnc * Sxyc - Sxc * Syc )/sqrt( (nnc*Sx2c - Sxc*Sxc) * (nnc*Sy2c - Syc*Syc));
summary.ccsum = Sxyc;

if(flag_par.verbose)
fprintf(fp1, "Pearson correlation coefficient: %6.3lf --- sum of products of charges %16.3lf e^2\n", summary.ccpc, summary.ccsum);
if(flag_par.verbose)
fclose(fp1);
summary.pppc = (nn * Sxy - Sx * Sy )/sqrt( (nn*Sx2 - Sx*Sx) * (nn*Sy2 - Sy*Sy));
summary.ppsum = Sxy;
if(flag_par.verbose && (flag_par.msms || flag_par.ns))
fprintf(fp2, "Pearson correlation coefficient: %6.3lf --- sum of products of pair potentials %16.3lf\n",
summary.pppc, summary.ppsum);
if(flag_par.msms || flag_par.ns)
for(m=0,i1=0,j1=-1; m< srf_ses_mol_1.n_srf_pt_tot; m++)
if(is_in_int_1[m])
{
if(i1%10000 == 0) j1++;
if(flag_par.verbose)
fprintf(fp3,"ATOM  %5i  SRF SRF %c%4i    %8.3lf%8.3lf%8.3lf  1.00%6.2lf      %s\n", 
i1%10000,
ch[j1],
i1%10000,
srf_ses_mol_1.srf_pt[m].r[0] , 
srf_ses_mol_1.srf_pt[m].r[1] , 
srf_ses_mol_1.srf_pt[m].r[2], srf_ses_mol_1.srf_pt[m].pot,
"    "
);
i1++;
}
if(flag_par.msms || flag_par.ns)
for(m=0,i2=0,j2=-1; m< srf_ses_mol_2.n_srf_pt_tot; m++)
if(is_in_int_2[m])
{
if(i2%10000 == 0) j2++;
if(flag_par.verbose)
fprintf(fp4,"ATOM  %5i  SRF SRF %c%4i    %8.3lf%8.3lf%8.3lf  1.00%6.2lf      %s\n", 
i2%10000,
ch[j2],
i2%10000,
srf_ses_mol_2.srf_pt[m].r[0] , 
srf_ses_mol_2.srf_pt[m].r[1] , 
srf_ses_mol_2.srf_pt[m].r[2], srf_ses_mol_2.srf_pt[m].pot,
"    "
);
i2++;
}

if(flag_par.verbose && (flag_par.msms || flag_par.ns))
fclose(fp3);
if(flag_par.verbose && (flag_par.msms || flag_par.ns))
fclose(fp4);
print_summary(summary,flag_par);
if(!flag_par.verbose)
{
if(flag_par.ns)
{
sprintf(buf_ns,"rm -f  %s_mol_1.prm %s_mol_2.prm %s.prm %s_mol_1.area %s_mol_2.area %s.area %s_mol_1.face %s_mol_2.face %s.face %s.vert %s_mol_1.vert %s_mol_2.vert  %s.xyzr %s_mol_1.xyzr %s_mol_2.xyzr",
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns,
flag_par.file_ns
);
if(system(buf_ns))
{ printf("the call: %s failed.... check it in your terminal\nexiting.....\n",buf_ns); exit(0);}
}
if(flag_par.msms)
{
sprintf(buf_ns,"rm -f  %s_mol_1.face %s_mol_2.face %s.face %s.vert %s_mol_1.vert %s_mol_2.vert  %s.xyzr %s_mol_1.xyzr %s_mol_2.xyzr %s.area %s.msms_log",
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms,
flag_par.file_msms
);
if(system(buf_ns))
{ printf("the call: %s failed.... check it in your terminal\nexiting.....\n",buf_ns); exit(0);}
}
}
}

void check_cmd_line(int argc, char *argv[], struct Flag_par *flag_par)
{
	int i;
        char tmp[100];
        char extension[100];

	if(argc < 5) 
	{
	printf("Usage:\n"); 
	printf("bluues2_cplx filename.pqr basename chain(s)_1 chain(s)_2 [Options]\n"); 
	printf("Options:\n"); 
	printf("-lite (output only summary file)\n"); 
	printf("-ns base_filename (base file name for NanoShaper) (base file name for NanoShaper parameter file)\n"); 
	printf("-rpns (read parameter file for NanoShaper instead of generating it)\n"); 
	printf("-m base_filename (base file name for msms)\n"); 
	printf("-pa x (area per point on SES in A^2, 0.1 default)\n"); 
	printf("-p x (use a probe radius for surfacing of x in A, 1.5 default)\n"); 
	printf("-tcr x (use a test charge radius of x in A, 0.7 default)\n"); 
	printf("-mr x (minimum atomic radius of x in A, 1.0 default. Lesser radii will be set to x A.)\n"); 
        printf("-s x (use a salt radius of x in A, 2.0 default)\n"); 
	printf("-pd x (inner dielectric constant, 20.0 default)\n"); 
	printf("-sd x (outer dielectric constant, 78.54 default)\n"); 
	printf("-kp x (factor for Still formula, 4 default)\n"); 
	printf("-i x (ionic strength (M), 0.15 default)\n"); 
	printf("-t x (temperature, 298.15 default)\n"); 
	printf("-c x (cutoff for GBR6 calculation of x A, 20.0 default\n"); 
	printf("-c2 x (cutoff short for averaging patches of x A, 8.0 default)\n"); 
	printf("-cs x (cutoff short for defining contacts between surfaces of x A, 1.0 default)\n"); 
	printf("-g x (surface tension coefficient kJ/(A^2 mol), 0.12 default)\n"); 
	printf("-srf (output surface points and potential)\n"); 
	printf("-srfpot (output atomic average surface potential (averaged at atomic SAS), surface and SES and SAS areas)\n"); 
	printf("-srfpotav (output atomic and average surface potential (averaged at atomic SAS over a neigborhood within c2 A, see above))\n"); 
	printf("\n"); 
	exit(1);
	}
        
        strcpy((*flag_par).file_in_pdb, argv[1]);
        strcpy((*flag_par).file_out, argv[2]);
        strcpy((*flag_par).chain_p,argv[3]);
        strcpy((*flag_par).chain_l,argv[4]);

	for (i = 5; i < argc; i++) 
        {
		if (!strcmp(argv[i],"-srfpot")) (*flag_par).print_srfpot = 1; 
		else if (!strcmp(argv[i],"-srf")) (*flag_par).print_surface = 1; 
		else if (!strcmp(argv[i],"-srfpotav")) (*flag_par).srfpotav = 1; 
		else if (!strcmp(argv[i],"-inf")) sscanf(argv[++i], "%lf", &((*flag_par).inflate)); 
		else if (!strcmp(argv[i],"-tcr")) sscanf(argv[++i], "%lf", &((*flag_par).test_charge_radius)); 
		else if (!strcmp(argv[i],"-pd")) sscanf(argv[++i], "%lf", &((*flag_par).pdie)); 
		else if (!strcmp(argv[i],"-sd")) sscanf(argv[++i], "%lf", &((*flag_par).sdie)); 
		else if (!strcmp(argv[i],"-mr")) sscanf(argv[++i], "%lf", &((*flag_par).min_radius)); 
                else if (!strcmp(argv[i],"-pa")) sscanf(argv[++i], "%lf", &((*flag_par).msms_area));
                else if (!strcmp(argv[i],"-dx")) (*flag_par).grid = 1;
		else if (!strcmp(argv[i],"-i")) sscanf(argv[++i], "%lf", &((*flag_par).ions)); 
		else if (!strcmp(argv[i],"-c2")) sscanf(argv[++i], "%lf", &((*flag_par).cutoff2)); 
		else if (!strcmp(argv[i],"-cs")) sscanf(argv[++i], "%lf", &((*flag_par).cutoff_short)); 
		else if (!strcmp(argv[i],"-kp")) sscanf(argv[++i], "%lf", &((*flag_par).kp)); 
 		else if (!strcmp(argv[i],"-m")) {sscanf(argv[++i], "%s", ((*flag_par).file_msms)); (*flag_par).msms = 1; (*flag_par).ns = 0;}
	        else if (!strcmp(argv[i],"-ns")) {sscanf(argv[++i], "%s", ((*flag_par).file_ns)); (*flag_par).ns = 1; (*flag_par).msms = 0;}
		else if (!strcmp(argv[i],"-rpns")) (*flag_par).rpns = 1; 
		else if (!strcmp(argv[i],"-p")) sscanf(argv[++i], "%lf", &((*flag_par).probe_radius)); 
		else if (!strcmp(argv[i],"-s")) sscanf(argv[++i], "%lf", &((*flag_par).salt_radius)); 
		else if (!strcmp(argv[i],"-t")) sscanf(argv[++i], "%lf", &((*flag_par).temp)); 
		else if (!strcmp(argv[i],"-g")) sscanf(argv[++i], "%lf", &((*flag_par).gamma)); 
		else if (!strcmp(argv[i],"-i")) sscanf(argv[++i], "%lf", &((*flag_par).ions)); 
		else if (!strcmp(argv[i],"-c")) sscanf(argv[++i], "%lf", &((*flag_par).cutoff)); 
		else if (!strcmp(argv[i],"-g")) (*flag_par).grid = 1; 
		else if (!strcmp(argv[i],"-lite")) (*flag_par).verbose = 0;
		else 
                {
                 printf("I don't know option %s\n", argv[i]);
                 exit(2);
                }
        }

}

void init_flag_par(struct Flag_par *flag_par)
{
(*flag_par).verbose=1;
(*flag_par).probe_radius=1.5;
(*flag_par).test_charge_radius=0.7;
(*flag_par).min_radius=1.0;
(*flag_par).inflate=0.0;
(*flag_par).salt_radius=2.0;
(*flag_par).ions=0.150;
(*flag_par).pdie=20.0;
(*flag_par).sdie=78.54;
(*flag_par).kp=4;
(*flag_par).temp=298.15;
(*flag_par).cutoff=20.0;
(*flag_par).cutoff2=8.0;
(*flag_par).cutoff_short=1.0;
(*flag_par).gbr = 1;
(*flag_par).srfpot = 1;
(*flag_par).print_srfpot = 0;
(*flag_par).print_surface = 0;
(*flag_par).srfpotav = 0;
(*flag_par).solv_e = 1;
(*flag_par).maxf = 999.94999;
(*flag_par).gamma = 0.12;
(*flag_par).grid = 0;
(*flag_par).ns=0;
strcpy((*flag_par).file_ns,"aux");
(*flag_par).rpns=0;
(*flag_par).msms=0;
(*flag_par).msms_area=0.5;
(*flag_par).pqg = 1;
(*flag_par).nx = 97;
(*flag_par).ny = 97;
(*flag_par).nz = 97;
(*flag_par).mesh = 1.0;
}

void print_info_flag_par(struct Flag_par flag_par)
{
        printf("--------------------------------\nRUN PARAMETERS:\n--------------------------------\n");
        printf("pdb file: %s\n", flag_par.file_in_pdb);
        printf("output file: %s\n", flag_par.file_out);
        printf("probe_radius: %lf\n", flag_par.probe_radius);
        printf("minimum atomic radius: %lf\n", flag_par.min_radius);
        printf("ionic strength (M): %lf\n", flag_par.ions);
        printf("salt_radius (A): %lf\n", flag_par.salt_radius);
        printf("inner dielectric constant: %lf\n", flag_par.pdie);
        printf("outer dielectric constant: %lf\n", flag_par.sdie);
        printf("temperature (K): %lf\n", flag_par.temp);
        printf("cutoff for GBR6 calculation (A): %lf\n", flag_par.cutoff);
        printf("cutoff for electrostatic patches (A): %lf\n", flag_par.cutoff2);
        printf("cutoff for contact definition between vdw surfs. (A): %lf\n", flag_par.cutoff_short);
        printf("gamma (kJ/A^2): %lf\n", flag_par.gamma);
        if(flag_par.srfpot || flag_par.grid)
        {
        printf("test charge radius: %lf\n", flag_par.test_charge_radius);
        }
        if(flag_par.msms )
        {
        printf("base file name for MSMS: %s\n", flag_par.file_msms);
        }
        if(flag_par.ns )
        {
        printf("base filename for NanoShaper: %s\n", flag_par.file_ns);
        }
        printf("area (A^2) per point on SES: %lf\n", flag_par.msms_area);
        printf("--------------------------------\n");
}

double fgb(double rgb1, double rgb2, double rw, double rs, double d, 
           double pdie, double sdie, double kd, double kp)
{
double f, rs_eff, rgb, deff;
rgb = sqrt(rgb1 * rgb2);
if(rs < rw) rs = rw; 
rs = sqrt(
(rgb1 + rs - rw) *
(rgb2 + rs - rw)
);
deff = sqrt(d*d + rgb * rgb *exp(-d*d/(kp * rgb * rgb)));


if(deff >= rs)
f = -1/(pdie * deff) + (1/sdie) * (exp(-kd * (deff - rs)) /(deff * (1 + kd * rs)));
else 
f = -1/(pdie * deff) + (1/sdie) * (1.0/(rs * (1 + kd * rs)) + 1.0/deff - 1.0/rs);

return f;
}

int cmp_int(const void *p1, const void *p2)
{
int i1, i2;
i1 = *((int *)p1);
i2 = *((int *)p2);
if(i2>i1) return -1;
else if(i1==i2) return 0;
else return 1;
}

int cmp_area_atom_by_string(const void *p1, const void *p2)
{
	struct Area_atom A, B;
	int check = 0 ; 

	A = *((struct Area_atom *)p1); 
	B = *((struct Area_atom *)p2); 

	if( A.n_model < B.n_model) check = -1; 
        else if ( strcmp(A.segment,B.segment) < 0) check = -1;
        else if ( A.chain < B.chain) check = -1;
        else if ( A.res_n < B.res_n) check = -1;
        else if ( A.res_ins < B.res_ins) check = -1;
        else if ( A.alt_loc < B.alt_loc) check = -1;
        else check = 1;
        return check;
}

int cmp_area_atom_by_n(const void *p1, const void *p2)
{
	struct Area_atom A, B;
	int check = 0 ; 

	A = *((struct Area_atom *)p1); 
	B = *((struct Area_atom *)p2); 

	if( A.at_n < B.at_n) check = -1; 
        else check = 1;
        return check;
}

void pqr2gbr6(struct System sys, struct Srf srf_ses, struct Srf srf_sas, struct Flag_par flag_par, struct Neighbour *neighbours, double *gbr6_ses, double *gbr6_sas)
{
int i, j, k, l;
double *phi,*ir1,*ir2,*ir3, dx = 1e-3;
double *v, tmpf, t0,t1,t2,t3, min[3], max[3], rave, *gbr6;
struct Srf *srf;

        v = (double *) calloc((size_t) 3,sizeof(double));
        phi=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir1=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir2=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
        ir3=(double *) calloc((size_t) sys.n_atoms, sizeof(double));
for(k=0; k<3; k++)
{
min[k] =  99999999999.9;
max[k] = -99999999999.9;
}


for(i=0; i<sys.n_atoms; i++)
for(k=0; k<3; k++)
{
if(sys.atoms[i].coor[k] < min[k]) min[k] = sys.atoms[i].coor[k];
if(sys.atoms[i].coor[k] > max[k]) max[k] = sys.atoms[i].coor[k];
}




rave = 0.5 * sqrt( 
       (max[0] - min[0])*(max[0] - min[0]) 
     + (max[1] - min[1])*(max[1] - min[1]) 
     + (max[2] - min[2])*(max[2] - min[2]) 
       )/sqrt(3.0);
if(rave < 1.0) rave = 1.0;

if(flag_par.msms || flag_par.ns)
	srf=&srf_ses;
else
	srf=&srf_sas;

for(i=0; i<sys.n_atoms; i++)
{
for(k=0; k <(*srf).at_srf[i].n_isrf; k++) 
{
l =(*srf).at_srf[i].isrf[k];
if(l >(*srf).n_srf_pt_tot) {printf("%i %i ....exiting....\n", l,(*srf).n_srf_pt_tot); exit(0);}

diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; 
if(tmpf > 0)
 {
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;

 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 


 }

}

for(j=0; j< neighbours[i].n_neighbours; j++) 
for(k=0; k <(*srf).at_srf[neighbours[i].list[j]].n_isrf; k++)
{
l =(*srf).at_srf[neighbours[i].list[j]].isrf[k];
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; 
 if(tmpf > 0)
 {
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;
 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
 }
}
}

for(i=0; i<sys.n_atoms; i++)
{
t0 = rave*rave*rave;
t1 = rave*t0;
t2 = rave*t1;
t3 = rave*t2;
if(phi[i] < 2*M_PI) phi[i] = 4*M_PI;
phi[i] = phi[i]/t0;
ir1[i] = ir1[i]/t1;
ir2[i] = ir2[i]/t2;
ir3[i] = ir3[i]/t3;

if(ir1[i] > 0.0)
{
ir1[i] = phi[i]/ir1[i];

}
else {
    ir1[i] = phi[i]/abs(ir1[i]);
      printf("Unexpected negative radius at atom  %5i %4s %4s %4i %c\n atom radius set to average radius, it may be due to approximate numerical treatment of internal cavities or surfaces\nThis will be treated approximately\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain);
     }
if(ir2[i] > 0) 
ir2[i] = pow(phi[i]/ir2[i], 1.0/2.0);

if(ir3[i] > 0) 
ir3[i] = pow(phi[i]/ir3[i], 1.0/3.0);

if(ir2[i] > ir1[i]) ir2[i] = ir1[i];
if(ir3[i] > ir2[i]) ir3[i] = ir2[i];
if(flag_par.msms || flag_par.ns)
	gbr6 = gbr6_ses;
	else
	gbr6 = gbr6_sas;

if(ir3[i] > 0) 
gbr6[i] = ir3[i];
else 
gbr6[i] = -1.0;
if(gbr6[i] < sys.atoms[i].radius) gbr6[i] = sys.atoms[i].radius;
if(gbr6[i] > flag_par.cutoff ) gbr6[i] = flag_par.cutoff;
}

if(!(flag_par.msms || flag_par.ns))
for(i=0; i<sys.n_atoms; i++)
{
if(gbr6[i] > 0.0)
{
if(gbr6[i] < sys.atoms[i].radius + flag_par.probe_radius) gbr6[i] = sys.atoms[i].radius +  flag_par.probe_radius;
if(gbr6[i] > flag_par.cutoff + flag_par.probe_radius) gbr6_ses[i] = flag_par.cutoff+ flag_par.probe_radius;
gbr6_ses[i] = (0.77322 + 0.06065 * flag_par.probe_radius) * gbr6[i] + (0.5142 -1.7585 * flag_par.probe_radius);
if(gbr6_ses[i] < sys.atoms[i].radius) gbr6_ses[i] = sys.atoms[i].radius;
if(gbr6_ses[i] > flag_par.cutoff ) gbr6_ses[i] = flag_par.cutoff;
}
else
{
gbr6[i] = DBL_MAX;
for(k=0; k< neighbours[i].n_neighbours; k++) 
{
j = neighbours[i].list[k];
if(gbr6[j] > 0.0 && gbr6[j] < gbr6[i]) 
gbr6[i] = gbr6[j];
}
}
if(gbr6[i] < sys.atoms[i].radius + flag_par.probe_radius) gbr6[i] = sys.atoms[i].radius +  flag_par.probe_radius;
if(gbr6[i] > flag_par.cutoff + flag_par.probe_radius) gbr6_ses[i] = flag_par.cutoff+ flag_par.probe_radius;
}
else
{

	srf = &srf_sas;
	gbr6 = gbr6_sas;
for(i=0; i<sys.n_atoms; i++)
{
 phi[i] = 0.0;
 ir1[i] = 0.0; 
 ir2[i] = 0.0; 
 ir3[i] = 0.0; 
for(k=0; k <(*srf).at_srf[i].n_isrf; k++) 
{
l =(*srf).at_srf[i].isrf[k];
if(l >(*srf).n_srf_pt_tot) {printf("%i %i ....exiting....\n", l,(*srf).n_srf_pt_tot); exit(0);}

diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);

tmpf = modv(v)/rave; 
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;



 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
}
for(j=0; j< neighbours[i].n_neighbours; j++) 
for(k=0; k <(*srf).at_srf[neighbours[i].list[j]].n_isrf; k++)
{
l =(*srf).at_srf[neighbours[i].list[j]].isrf[k];
diffv(v,(*srf).srf_pt[l].r, sys.atoms[i].coor);
tmpf = modv(v)/rave; 
 t0 =(*srf).srf_pt[l].a*dotv(v,(*srf).srf_pt[l].vec)/(tmpf*tmpf*tmpf);
 t1 = t0/tmpf;
 t2 = t1/tmpf;
 t3 = t2/tmpf;
 phi[i] = phi[i] + t0;
 ir1[i] = ir1[i] + t1; 
 ir2[i] = ir2[i] + t2; 
 ir3[i] = ir3[i] + t3; 
}


}
for(i=0; i<sys.n_atoms; i++)
{
t0 = rave*rave*rave;
t1 = rave*t0;
t2 = rave*t1;
t3 = rave*t2;
phi[i] = phi[i]/t0;
ir1[i] = ir1[i]/t1;
ir2[i] = ir2[i]/t2;
ir3[i] = ir3[i]/t3;

if(ir1[i] > 0.0)
{
ir1[i] = phi[i]/ir1[i];

}
else {
      ir1[i] = rave;
      printf("Unexpected negative radius at atom  %5i %4s %4s %4i %c \n atom radius set to average radius, it may be due to approximate numerical treatment of internal cavities\nTry with msms surface\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain);
     }
if(ir2[i] > 0) 
ir2[i] = pow(phi[i]/ir2[i], 1.0/2.0);
else ir2[i] = ir1[i] + dx;
if(ir3[i] > 0) 
ir3[i] = pow(phi[i]/ir3[i], 1.0/3.0);
else ir3[i] = ir2[i] + dx;
if(ir2[i] > ir1[i]) ir2[i] = ir1[i];
if(ir3[i] > ir2[i]) ir3[i] = ir2[i];
gbr6[i] = ir3[i];
}
}
}

void srfpot(struct System *sys, struct Srf *srf_ses, struct Srf *srf_sas, double *gbr6, struct Neighbour *neighbours, struct Flag_par flag_par)
{
int i, j, k, l, jj;
double d, f, kd, k_el;
double      q = 1.602176487E-19;
double      N_av = 6.02214179E23;
double      e0 = 8.8541878176E-12;
double      kb = 1.3806504E-23;
double      RT = 1.3806504 * 6.02214179 * flag_par.temp;
struct Srf *srf;
if(flag_par.msms || flag_par.ns) srf = srf_ses; else srf = srf_sas;
k_el = N_av * 1E10 * q * q /(4 * M_PI * e0);
kd = 1E-10 * sqrt((1000 * 2 * flag_par.ions * N_av * q * q )/
                         (flag_par.sdie * e0 * kb * flag_par.temp));
for(l=0; l< (*srf).n_srf_pt_tot; l++) (*srf).srf_pt[l].pot = 0.0;
for(i=0; i<(*sys).n_atoms; i++)
{
(*sys).atoms[i].temp = 0.0;
for(k=0; k < (*srf).at_srf[i].n_isrf; k++) 
{
l = (*srf).at_srf[i].isrf[k];
d = distv((*srf).srf_pt[l].r, (*sys).atoms[i].coor);
f = fgb(gbr6[i],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  kd, flag_par.kp);
(*srf).srf_pt[l].pot = (*srf).srf_pt[l].pot + (*sys).atoms[i].charge/(d *flag_par.pdie) + f * (*sys).atoms[i].charge;
(*sys).atoms[i].temp = (*sys).atoms[i].temp + (*sys).atoms[i].charge/(d *flag_par.pdie); 
(*sys).atoms[i].temp = (*sys).atoms[i].temp + f * (*sys).atoms[i].charge; 
for(j=0; j< neighbours[i].n_neighbours; j++) 
{
jj = neighbours[i].list[j];
d = distv((*srf).srf_pt[l].r, (*sys).atoms[jj].coor);
f = fgb(gbr6[jj],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  kd, flag_par.kp);
(*srf).srf_pt[l].pot = (*srf).srf_pt[l].pot + (*sys).atoms[jj].charge/(d *flag_par.pdie) + f * (*sys).atoms[jj].charge;
(*sys).atoms[i].temp = (*sys).atoms[i].temp + (*sys).atoms[jj].charge/(d *flag_par.pdie); 
(*sys).atoms[i].temp = (*sys).atoms[i].temp + f * (*sys).atoms[jj].charge; 
}
}
if((*srf).at_srf[i].n_isrf != 0)
{
(*sys).atoms[i].temp = (k_el/1E3) * (*sys).atoms[i].temp / (double) (*srf).at_srf[i].n_isrf;
(*sys).atoms[i].occ = 1.0;
}
else 
{
(*sys).atoms[i].temp = 0.0;
(*sys).atoms[i].occ = 0.0;
}
}

for(l=0; l< (*srf).n_srf_pt_tot; l++) (*srf).srf_pt[l].pot =  (k_el/1E3) * (*srf).srf_pt[l].pot;
}

void scalev(double t, double *v)
{
    v[0] = t*v[0];
    v[1] = t*v[1];
    v[2] = t*v[2];
}

void get_srf(struct System sys, struct Atom_grid *patom_grid, struct Srf *psrf_ses, struct Srf *psrf_sas, struct Flag_par flag_par)
{
      char buf[2048];
      FILE *fp;
      int i,j,k,l,n_faces, max_n_srf_pt;
      struct Srf tmp_srf, srf;
      struct Atom_grid atom_grid;
      double d1, d2, d3, sp, At, A, t, r_max;
      struct Msms_face *faces;
      
      atom_grid = *patom_grid;

        for(i=0; i< sys.n_atoms; i++)
        if((sys.atoms[i].radius + flag_par.probe_radius) > r_max)
            r_max = sys.atoms[i].radius + flag_par.probe_radius;
       (*psrf_ses).at_srf = (struct At_srf *) calloc( (size_t) sys.n_atoms, sizeof(struct At_srf));
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_ses).at_srf[i].n_isrf=0;

        (*psrf_sas).at_srf = (struct At_srf *) calloc( (size_t) sys.n_atoms , sizeof(struct At_srf));
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_sas).at_srf[i].n_isrf=0;
if(flag_par.msms || flag_par.ns)
{
if(flag_par.msms) 
{
printf("Calling msms for molecular surface calculation...\n");
printf("Make sure that msms is in your PATH...\n");
printf("and that the minimum radius is set as to avoid one atom included in one another (e.g. use option -mr 1.0)\n");
sprintf(buf,"%s.xyzr",flag_par.file_msms);
fp = file_open(buf,"w");
for(i = 0; i< sys.n_atoms; i++)
fprintf(fp,"%11.3f %11.3f %11.3f %11.3f\n", 
sys.atoms[i].coor[0],
sys.atoms[i].coor[1],
sys.atoms[i].coor[2],
sys.atoms[i].radius);
fclose(fp);
sprintf(buf,"msms -if %s.xyzr -of %s -af %s -density %lf -probe_radius %lf -no_rest > %s.msms_log\n",flag_par.file_msms, flag_par.file_msms, flag_par.file_msms, 1.0/flag_par.msms_area, flag_par.probe_radius, flag_par.file_msms);
printf("running... %s\n", buf);
if(system(buf)) { printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n",buf); exit(0);}
sprintf(buf,"%s.vert",flag_par.file_msms);
fp = file_open(buf,"r");
printf(".....I am reading surface points and normal vectors in %s.vert\n", flag_par.file_msms);
}
else if(flag_par.ns)
{
sprintf(buf,"%s.xyzr",flag_par.file_ns);
fp = file_open(buf,"w");
for(i = 0; i< sys.n_atoms; i++)
fprintf(fp,"%11.3f %11.3f %11.3f %11.3f\n", 
sys.atoms[i].coor[0],
sys.atoms[i].coor[1],
sys.atoms[i].coor[2],
sys.atoms[i].radius);
fclose(fp);


sprintf(buf,"%s.prm",flag_par.file_ns);
if(!flag_par.rpns)
{
fp = file_open(buf,"w");
fprintf(fp,"Operative_Mode = normal\n");
fprintf(fp,"Compute_Vertex_Normals = true\n");
fprintf(fp,"Save_Mesh_MSMS_Format = true\n");
fprintf(fp,"Grid_scale = %lf\n", 1/sqrt(1.6 * flag_par.msms_area));
fprintf(fp,"Grid_perfil = %lf\n", (90 * (double) sys.n_atoms / (double) (90 + sys.n_atoms) ) );
fprintf(fp,"XYZR_FileName = %s.xyzr\n", flag_par.file_ns);
fprintf(fp,"Build_epsilon_maps  = false\n");
fprintf(fp,"Build_status_map = true\n");
fprintf(fp,"Surface = ses\n");
fprintf(fp,"Smooth_Mesh = false\n");
fprintf(fp,"Number_thread = -1 #default value is 32\n");
fprintf(fp,"Probe_Radius = %lf\n", flag_par.probe_radius);
fprintf(fp,"Self_Intersections_Grid_Coefficient = %lf\n", flag_par.probe_radius + 1.0);
fprintf(fp,"Accurate_Triangulation = true\n");
fprintf(fp,"Triangulation = true\n");
fprintf(fp,"Check_duplicated_vertices = true\n");
fprintf(fp,"Vertex_Atom_Info = true\n");
fclose(fp);
}
//sprintf(buf,"$(ls -l $(which NanoShaper) | awk '{print $NF}') %s.prm\ncp triangulatedSurf.face %s.face\ncp triangulatedSurf.vert %s.vert\ncp triangleAreas.txt %s.area\n",flag_par.file_ns,flag_par.file_ns,flag_par.file_ns,flag_par.file_ns);
sprintf(buf,"NanoShaper %s.prm\ncp triangulatedSurf.face %s.face\ncp triangulatedSurf.vert %s.vert\ncp triangleAreas.txt %s.area\n",flag_par.file_ns,flag_par.file_ns,flag_par.file_ns,flag_par.file_ns);
printf("running NanoShaper... %s\n", buf);
if(system(buf)) { printf("the call: %s did not succeed, check it in your terminal\nexiting.....\n",buf); exit(0);}
sprintf(buf,"%s.vert",flag_par.file_ns);
fp = file_open(buf,"r");
printf(".....I am reading surface points and normal vectors in %s.vert\n", flag_par.file_ns);
}
do 
{
if(fgets(buf,120,fp) != NULL) ;
else 
{printf("fgets(buf,120,fp) returned NULL...., exiting\n"); exit(0);}
} 
while (buf[0] == '#');
sscanf(buf,"%i", &(tmp_srf.n_srf_pt_tot));
tmp_srf.srf_pt= (struct Srf_pt *) calloc((size_t) tmp_srf.n_srf_pt_tot, sizeof(struct Srf_pt));
        for(i=0; i < tmp_srf.n_srf_pt_tot ; i++)
        {
        tmp_srf.srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
        tmp_srf.srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
        }

k=0;
while(fgets(buf, 120, fp) != NULL)
{

sscanf(buf,"%lf %lf %lf %lf %lf %lf %i %i %i\n",
&(tmp_srf.srf_pt[k].r[0]),
&(tmp_srf.srf_pt[k].r[1]),
&(tmp_srf.srf_pt[k].r[2]),
&(tmp_srf.srf_pt[k].vec[0]),
&(tmp_srf.srf_pt[k].vec[1]),
&(tmp_srf.srf_pt[k].vec[2]),
&(tmp_srf.srf_pt[k].atn), 
&j,
&(tmp_srf.srf_pt[k].type));
tmp_srf.srf_pt[k].at = j - 1;
tmp_srf.srf_pt[k].a= flag_par.msms_area;

k++;
if(k > tmp_srf.n_srf_pt_tot) {printf("more points in file than stated in header...\nexiting...\n"); exit(0);}
}
fclose(fp);
tmp_srf.n_srf_pt_tot = k;
if(flag_par.msms)
{
printf("I have read %i surface points in %s.vert\n", tmp_srf.n_srf_pt_tot, flag_par.file_msms);
sprintf(buf,"%s.face",flag_par.file_msms);
fp = file_open(buf,"r");
printf(".....I am reading surface faces in %s.face\n", flag_par.file_msms);
}
else if(flag_par.ns)
{
printf("I have read %i surface points in %s.vert\n", tmp_srf.n_srf_pt_tot, flag_par.file_ns);
sprintf(buf,"%s.face",flag_par.file_ns);
fp = file_open(buf,"r");
printf(".....I am reading surface faces in %s.face\n", flag_par.file_ns);
}
do 
{
if(fgets(buf,120,fp) != NULL) ;
else 
{printf("fgets(buf,120,fp) returned NULL...., exiting\n"); exit(0);}
} 
while (buf[0] == '#');
sscanf(buf,"%i", &n_faces);
faces = (struct Msms_face *) calloc((size_t) n_faces, sizeof(struct Msms_face));
k=0;
while(fgets(buf, 120, fp) != NULL)
{
sscanf(buf,"%i %i %i %i %*i\n",
&(faces[k].n1),
&(faces[k].n2),
&(faces[k].n3),
&(faces[k].type)
);
faces[k].n1 = faces[k].n1 - 1;
faces[k].n2 = faces[k].n2 - 1;
faces[k].n3 = faces[k].n3 - 1;
k++;
if(k > n_faces) {printf("more points in file than stated in header...\nexiting...\n"); exit(0);}
}
n_faces = k;
if(flag_par.msms)
printf("I have read %i surface faces in %s.face\n", n_faces, flag_par.file_msms);
else if(flag_par.ns)
printf("I have read %i surface faces in %s.face\n", n_faces, flag_par.file_ns);
j = 0;
for(k = 0; k< n_faces; k++)
{

if((tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at)  &&
(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at) )
{
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n3].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else if(tmp_srf.srf_pt[faces[k].n2].at == tmp_srf.srf_pt[faces[k].n3].at)
{ 
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
else
{
l = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
l = tmp_srf.srf_pt[faces[k].n3].at;
(*psrf_ses).at_srf[l].n_isrf++;
j++;
}
}
(*psrf_ses).n_srf_pt_tot = j;
        (*psrf_ses).srf_pt= (struct Srf_pt*) calloc((size_t) (*psrf_ses).n_srf_pt_tot, sizeof(struct Srf_pt));
        for(i=0; i < (*psrf_ses).n_srf_pt_tot ; i++)
         {
         (*psrf_ses).srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
         (*psrf_ses).srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
         }



j = 0;
for(k = 0, At = 0.0; k< n_faces; k++)
{
for(i = 0; i < 3 ; i++)
{
(*psrf_ses).srf_pt[j].r[i] = 
 (tmp_srf.srf_pt[faces[k].n1].r[i] + 
  tmp_srf.srf_pt[faces[k].n2].r[i] + 
  tmp_srf.srf_pt[faces[k].n3].r[i] )/3.0; 
(*psrf_ses).srf_pt[j].vec[i] = 
 (tmp_srf.srf_pt[faces[k].n1].vec[i] + 
  tmp_srf.srf_pt[faces[k].n2].vec[i] + 
  tmp_srf.srf_pt[faces[k].n3].vec[i] ); 
}
t = modv((*psrf_ses).srf_pt[j].vec);
for(i = 0; i < 3 ; i++)
  (*psrf_ses).srf_pt[j].vec[i]  = (*psrf_ses).srf_pt[j].vec[i] / t; 
d1 = distv(tmp_srf.srf_pt[faces[k].n1].r, tmp_srf.srf_pt[faces[k].n2].r);
d2 = distv(tmp_srf.srf_pt[faces[k].n2].r, tmp_srf.srf_pt[faces[k].n3].r);
d3 = distv(tmp_srf.srf_pt[faces[k].n3].r, tmp_srf.srf_pt[faces[k].n1].r);
sp = (d1 + d2 + d3)/2.0;
A = sqrt(sp*(sp-d1)*(sp-d2)*(sp-d3));
  (*psrf_ses).srf_pt[j].a  = A; 
At = At + A;


(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
if(tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at  &&
tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at ) 
{
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).srf_pt[j].a = A;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3 )
(*psrf_ses).srf_pt[j].type = -1;
(*psrf_ses).srf_pt[j].type = -2;
j++;
}

else if(
   tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n2].at != tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n1].at;
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n3].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type ;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}

else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n1].at == tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}

else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n2].at == tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].a = 2*A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}

else if(
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n2].at &&
   tmp_srf.srf_pt[faces[k].n1].at != tmp_srf.srf_pt[faces[k].n3].at &&
   tmp_srf.srf_pt[faces[k].n2].at != tmp_srf.srf_pt[faces[k].n3].at) 
{
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n2].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
(*psrf_ses).srf_pt[j].a = A/3;
(*psrf_ses).srf_pt[j].at = tmp_srf.srf_pt[faces[k].n3].at;
for(i = 0; i < 3 ; i++)
  {
  (*psrf_ses).srf_pt[j].r[i] = (*psrf_ses).srf_pt[j-1].r[i];
  (*psrf_ses).srf_pt[j].vec[i] = (*psrf_ses).srf_pt[j-1].vec[i];
  }
(*psrf_ses).srf_pt[j].type = faces[k].type;
if(faces[k].type != 3)
(*psrf_ses).srf_pt[j].type = -1;
j++;
}
}

for(i = 0; i <  sys.n_atoms ; i++)
{
(*psrf_ses).at_srf[i].isrf = (int *) calloc((size_t) ((*psrf_ses).at_srf[i].n_isrf), sizeof(int));
(*psrf_ses).at_srf[i].n_isrf = 0;
}

for(i = 0; i <  (*psrf_ses).n_srf_pt_tot ; i++)
 (*psrf_ses).at_srf[(*psrf_ses).srf_pt[i].at].isrf[(*psrf_ses).at_srf[(*psrf_ses).srf_pt[i].at].n_isrf++] = i;
        for(i=0, A=0; i < (*psrf_ses).n_srf_pt_tot ; i++)
        A = (*psrf_ses).srf_pt[i].a + A;
        (*psrf_ses).area = A;
       printf("Total SES area: %lf A^2\n", A);
}

{

max_n_srf_pt = (int) (sys.n_atoms * 4.0 * M_PI * r_max * r_max/ flag_par.msms_area );


        (*psrf_sas).srf_pt= (struct Srf_pt *) calloc((size_t) max_n_srf_pt, sizeof(struct Srf_pt));

        for(i=0; i < max_n_srf_pt ; i++)
         {
         (*psrf_sas).srf_pt[i].r=(double *) calloc((size_t) 3, sizeof(double));
         (*psrf_sas).srf_pt[i].vec=(double *) calloc((size_t) 3, sizeof(double));
         }
        for(i=0; i< sys.n_atoms; i++)
        (*psrf_sas).at_srf[i].isrf= (int *) calloc( (size_t) (max_n_srf_pt/ sys.n_atoms) , sizeof(int));
        (*psrf_sas).alloc_srf_pt = max_n_srf_pt;


grid_srf_dens(sys, atom_grid, flag_par, psrf_sas);
}

}
void print_radii(struct System sys, struct Trj trj, double *gbr6_ses, double *gbr6_sas, struct Flag_par flag_par)
{
	double *tmpr;
	int i;
	char *buf;
	FILE *fp;
buf= (char *) calloc((size_t) 1024, sizeof(char));
if(flag_par.pqg)
{
tmpr = (double *) calloc((size_t) sys.n_atoms, sizeof(double));
for(i=0; i<sys.n_atoms; i++)
	tmpr[i] = sys.atoms[i].radius;
for(i=0; i<sys.n_atoms; i++)
{

if(gbr6_ses[i] < sys.atoms[i].radius)
	gbr6_ses[i] = sys.atoms[i].radius;
sys.atoms[i].occ = sys.atoms[i].charge;
sys.atoms[i].temp = gbr6_ses[i];
}
for(i = 0; i < sys.n_atoms; i++)
tmpr[i] = sys.atoms[i].radius;
if(flag_par.verbose)
{
strcpy(buf,flag_par.file_out);
strcat(buf,".pqg");
fp = fopen(buf,"w");
printf(".....I am writing pqg file in\n%s\n",buf);
write_PDB_atoms(fp, sys.n_atoms, sys.atoms, trj);
fclose(fp);
}

}
for(i=0; i<sys.n_atoms; i++)
sys.atoms[i].radius = tmpr[i]; 
for(i=0; i<sys.n_atoms; i++)
{
if(flag_par.msms || flag_par.ns) 
{
if( gbr6_ses[i] < sys.atoms[i].radius)
gbr6_ses[i] = sys.atoms[i].radius;
}

else if( gbr6_sas[i] < (sys.atoms[i].radius + flag_par.probe_radius))
gbr6_sas[i] = sys.atoms[i].radius + flag_par.probe_radius;



}

strcpy(buf,flag_par.file_out);
strcat(buf,".gbr");
fp = fopen(buf,"w");
        printf(".....I am writing GB radii in:\n%s\n",buf);
fprintf(fp, "       N  ATOM  RES RESN CH.   GBR6(SES)  GBR6(SAS)\n");
for(i=0; i<sys.n_atoms; i++)
   fprintf(fp, "GBR %5i %4s %4s %4i %c  %9.3f  %9.3f\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               gbr6_ses[i],
               gbr6_sas[i]
               );
fclose(fp);
        printf("--------------------------------\n");
free(tmpr);
free(buf);
}





void print_srf_pot(struct System sys, struct Trj trj, double *gbr6, struct Srf *srf_ses, struct Srf *srf_sas, struct Neighbour *neighbours, struct Area_atom *area_atom, struct Flag_par flag_par)
{
      char buf[1024];
      FILE *fp;
      int i,j,jj;
      struct System model;
      double tmpf;
      char chain_old, res_ins_old, res_name_old[5];
      int res_n_old;
      double area, area_ses,area_sas;
      struct Srf *srf;
char *ch="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
char segid[5];
 
srf = srf_sas;
if(flag_par.srfpot )
{    
if(flag_par.print_srfpot && flag_par.verbose)
{
strcpy(buf,flag_par.file_out);
strcat(buf,".srfatpot");
fp = fopen(buf,"w");
printf(".....I am writing atomic surface potentials in\n%s\n",buf);
write_PDB_atoms(fp, sys.n_atoms, sys.atoms, trj);
fclose(fp);
        printf("--------------------------------\n");
}

model.atoms = (struct Atom *) calloc((size_t) sys.n_atoms, sizeof(struct Atom));
for(i=0; i< sys.n_atoms; i++)
area_atom[i].area = 0.0; 
for(i=0; i<(*srf).n_srf_pt_tot; i++)
area_atom[(*srf).srf_pt[i].at].area = area_atom[(*srf).srf_pt[i].at].area + (*srf).srf_pt[i].a;
for(i=0; i< sys.n_atoms; i++)
copy_atom(&(model.atoms[i]),sys.atoms[i]);
model.n_atoms = sys.n_atoms;

if(flag_par.srfpotav)
{
for(i=0; i< model.n_atoms; i++)
{
tmpf=area_atom[i].area; 
model.atoms[i].temp = sys.atoms[i].temp * area_atom[i].area;
for(jj=0; jj<neighbours[i].n_neighbours; jj++)
{
j = neighbours[i].list[jj];
if(distv(model.atoms[i].coor,model.atoms[j].coor) <=  flag_par.cutoff2 && i != j)
{
tmpf = tmpf + area_atom[j].area;
model.atoms[i].temp = model.atoms[i].temp + area_atom[j].area * model.atoms[j].temp;
}
}
if(area_atom[i].area > 0.0)
model.atoms[i].temp = model.atoms[i].temp/tmpf;
else
model.atoms[i].temp = 0.0;
}
strcpy(buf,flag_par.file_out);
strcat(buf,".srfatpotav");
fp = fopen(buf,"w");
printf(".....I am writing averaged atomic surface potentials in\n%s\n",buf);
write_PDB_atoms(fp, model.n_atoms, model.atoms, trj);
fclose(fp);
}
//        printf("--------------------------------\n");
free(model.atoms);
}
      
if(flag_par.print_surface)
{
	if(flag_par.msms || flag_par.ns)
srf = srf_ses;
	else
srf = srf_sas;


if(flag_par.verbose)
{
printf(".....I am writing %i surface points in\n%s\n",(*srf).n_srf_pt_tot,buf);
strcpy(buf,flag_par.file_out);
strcat(buf,".at_r_n");
fp = fopen(buf,"w");
printf(".....I am writing surface point atom n., coordinate and normal versor in\n%s\n",buf);
fprintf(fp,"ATNUM     X        Y        Z        VX       VY       VZ   \n");
for(i=0; i<(*srf).n_srf_pt_tot; i++)
fprintf(fp,"%5i %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
sys.atoms[(*srf).srf_pt[i].at].at_n,
(*srf).srf_pt[i].r[0] , 
(*srf).srf_pt[i].r[1] , 
(*srf).srf_pt[i].r[2] ,
(*srf).srf_pt[i].vec[0] , 
(*srf).srf_pt[i].vec[1] , 
(*srf).srf_pt[i].vec[2] 
);
fclose(fp);
        printf("--------------------------------\n");
}

if(flag_par.verbose)
{
strcpy(buf,flag_par.file_out);
strcat(buf,".area");
printf(".....I am writing surface area by sys, chain, residue, atom in\n%s\n",buf);
fp = fopen(buf,"w");
}
for(i=0; i<sys.n_atoms; i++)
area_atom[i].area = area_atom[i].ses = area_atom[i].sas = 0.0;
if(flag_par.msms || flag_par.ns)
{
srf = srf_ses;
for(i=0; i<(*srf).n_srf_pt_tot; i++)
area_atom[(*srf).srf_pt[i].at].ses = area_atom[(*srf).srf_pt[i].at].ses + (*srf).srf_pt[i].a; 
}
srf = srf_sas;
for(i=0; i<(*srf).n_srf_pt_tot; i++)
area_atom[(*srf).srf_pt[i].at].sas = area_atom[(*srf).srf_pt[i].at].sas + (*srf).srf_pt[i].a; 
for(i=0; i<sys.n_atoms; i++)
{
strcpy(area_atom[i].at_name,sys.atoms[i].at_name);
area_atom[i].alt_loc = sys.atoms[i].alt_loc;
strcpy(area_atom[i].res_name, sys.atoms[i].res_name);
area_atom[i].chain = sys.atoms[i].chain;
if(area_atom[i].chain == ' ') area_atom[i].chain = '_';
area_atom[i].at_n = i+1;
area_atom[i].res_n = sys.atoms[i].res_n;
area_atom[i].res_ins = sys.atoms[i].res_ins;
}

qsort(area_atom, sys.n_atoms, sizeof(struct Area_atom), &cmp_area_atom_by_string);
for(i=0, area_ses = area_sas = 0.0; i<sys.n_atoms; i++)
{
area_ses = area_ses + area_atom[i].ses;
area_sas = area_sas + area_atom[i].sas;
}


if(flag_par.verbose)
fprintf(fp, "System area = %.5e A^2 (SES)  %.5e A^2 (SAS)\n", area_ses, area_sas);
(*srf_ses).area = area_ses;
(*srf_sas).area = area_sas;

for(i=0,area_ses = area_sas = 0.0; i<sys.n_atoms; i++)
{
if(i == 0 || area_atom[i].chain == chain_old)
{
area_ses = area_ses + area_atom[i].ses;
area_sas = area_sas + area_atom[i].sas;
}
else 
{
if(flag_par.verbose)
fprintf(fp, "Chain %c area = %.5e A^2 (SES)  %.5e A^2 (SAS)\n", chain_old, area_ses, area_sas);
area_ses = area_atom[i].ses;
area_sas = area_atom[i].sas;
}
chain_old = area_atom[i].chain;
if(i == sys.n_atoms - 1)
if(flag_par.verbose)
fprintf(fp, "Chain %c area = %.5e A^2 (SES)  %.5e A^2 (SAS)\n", chain_old, area_ses, area_sas);
}


for(i=0,area_ses = area_sas = 0.0; i<sys.n_atoms; i++)
{
if(i == 0 || (area_atom[i].chain == chain_old && 
area_atom[i].res_n == res_n_old && area_atom[i].res_ins == res_ins_old))
{
area_ses = area_ses + area_atom[i].ses;
area_sas = area_sas + area_atom[i].sas;
}
else 
{
if(flag_par.verbose)
fprintf(fp, "Residue %4s %5i %c chain %c area = %.4e A^2 (SES)   %.4e A^2 (SAS)\n", res_name_old, res_n_old, res_ins_old, chain_old, area_ses, area_sas);
area_ses = area_atom[i].ses;
area_sas = area_atom[i].sas;
}
chain_old = area_atom[i].chain;
res_n_old = area_atom[i].res_n;
res_ins_old = area_atom[i].res_ins;
strcpy(res_name_old, area_atom[i].res_name);
if(i == sys.n_atoms - 1)
if(flag_par.verbose)
fprintf(fp, "Residue %4s %5i %c chain %c area = %.4e A^2 (SES)   %.4e A^2 (SAS)\n", res_name_old, res_n_old, res_ins_old, chain_old, area_ses, area_sas);
}
qsort(area_atom, sys.n_atoms, sizeof(struct Area_atom), &cmp_area_atom_by_n);
if(flag_par.verbose)
for(i=0; i<sys.n_atoms; i++)
fprintf(fp, "Atom %4s %4s %5i %c chain %c area = %.3e A^2 (SES)  %.3e A^2 (SAS)\n", area_atom[i].at_name, area_atom[i].res_name, area_atom[i].res_n, area_atom[i].res_ins, area_atom[i].chain, area_atom[i].ses, area_atom[i].sas);

if(flag_par.verbose)
fclose(fp);

if(flag_par.verbose)
        printf("--------------------------------\n");
}
}


void calc_nrg(struct System sys, struct Energy *pnrg, double *gbr6, struct Srf srf, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par)
{
int i, j, jj;
double t1, d, *v, f, *tmpgbr;
struct Energy nrg;
char buf[120];
FILE *fp;
v=(double *) calloc((size_t) 3, sizeof(double));

if(flag_par.solv_e || 1)  
{
(*pnrg).area=0; 
(*pnrg).pol=0; 
(*pnrg).born = 0;
(*pnrg).coul=0; 
(*pnrg).solv=0; 
(*pnrg).tot=0; 
for(i=0; i<(srf).n_srf_pt_tot; i++)
{
(*pnrg).area = (*pnrg).area + (srf).srf_pt[i].a;
(*pnrg).area_a[(srf).srf_pt[i].at] = (*pnrg).area_a[(srf).srf_pt[i].at] + (srf).srf_pt[i].a;
}
for(i=0; i<sys.n_atoms; i++)
(*pnrg).area_a[i] = flag_par.gamma * (*pnrg).area_a[i]; 
(*pnrg).area = (*pnrg).area * flag_par.gamma;

for(i=0; i<sys.n_atoms; i++)
{
for(jj=0; jj<neighbours[i].n_neighbours; jj++)
if(neighbours[i].list[jj] > i)
{
j = neighbours[i].list[jj];
t1 = constants.k_el * sys.atoms[i].charge * sys.atoms[j].charge/1000;
diffv(v,sys.atoms[i].coor, sys.atoms[j].coor);
d = modv(v);
f = fgb(gbr6[i],gbr6[j],  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
(*pnrg).pol_a[i] = (*pnrg).pol_a[i] + t1 * f/2.0;
(*pnrg).pol_a[j] = (*pnrg).pol_a[j] + t1 * f/2.0;
(*pnrg).pol = (*pnrg).pol + t1 * f;
if(d<=0) {printf("Distance between atoms %i and %i < 0, exiting....\n", i, j); 
printf("coord %i: %lf %lf %lf\ncoord %i: %lf %lf %lf\n", 
i,
sys.atoms[i].coor[0],
sys.atoms[i].coor[1],
sys.atoms[i].coor[2],
j,
sys.atoms[j].coor[0],
sys.atoms[j].coor[1],
sys.atoms[j].coor[2]);
 exit(0);}
(*pnrg).coul_a[i] = (*pnrg).coul_a[i] + t1/(2.0 * flag_par.pdie * d);
(*pnrg).coul_a[j] = (*pnrg).coul_a[j] + t1/(2.0 * flag_par.pdie * d);
(*pnrg).coul = (*pnrg).coul + t1/(flag_par.pdie * d);
}

d = 0;
t1 = constants.k_el * sys.atoms[i].charge * sys.atoms[i].charge/1000.0;

f = fgb(gbr6[i],gbr6[i],  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
(*pnrg).born = (*pnrg).born + t1 * f / 2.0;

(*pnrg).born_a[i] = (*pnrg).born_a[i] + t1 * f / 2.0 ;
}

}
//if(!flag_par.msms && !flag_par.ns)
//for(i=0; i<sys.n_atoms; i++)
//gbr6[i] = tmpgbr[i];
for(i=0; i<sys.n_atoms; i++)
{
(*pnrg).solv_a[i] = (*pnrg).born_a[i] + (*pnrg).pol_a[i] + (*pnrg).area_a[i];
(*pnrg).tot_a[i] = (*pnrg).coul_a[i] + (*pnrg).solv_a[i];
}
(*pnrg).solv = (*pnrg).born + (*pnrg).pol + (*pnrg).area;
(*pnrg).tot = (*pnrg).coul + (*pnrg).solv;
}

void print_nrg(struct Energy nrg, struct System sys, struct Flag_par flag_par)
{
	FILE *fp;
	int i;
	char buf[1024];
printf("ENERGIES\n");
printf("--------------------------------\n");
printf("Total Coulomb   energy: %14.4lf (kJ/mol)\n", nrg.coul);
printf("Total solvation energy: %14.4lf (kJ/mol)\n", nrg.solv);
printf(" --Born self    energy: %14.4lf (kJ/mol)\n", nrg.born);
printf(" --Coulomb solv energy: %14.4lf (kJ/mol)\n", nrg.pol);
printf(" --sasa         energy: %14.4lf (kJ/mol)\n", nrg.area);
printf("Total           energy: %14.4lf (kJ/mol)\n", (nrg.coul + nrg.solv));
printf("--------------------------------\n");
 
strcpy(buf,flag_par.file_out);
strcat(buf,".nrg");
fp = fopen(buf,"w");
printf(".....I am writing total and atomic energies in\n%s\n",buf);
fprintf(fp,"Total Coulomb   energy: %14.4lf (kJ/mol)\n", nrg.coul);
fprintf(fp,"Total solvation energy: %14.4lf (kJ/mol)\n", nrg.solv);
fprintf(fp," --Born self    energy: %14.4lf (kJ/mol)\n", nrg.born);
fprintf(fp," --Coulomb solv energy: %14.4lf (kJ/mol)\n", nrg.pol);
fprintf(fp," --sasa         energy: %14.4lf (kJ/mol)\n", nrg.area);
fprintf(fp,"Total           energy: %14.4lf (kJ/mol)\n", (nrg.coul + nrg.solv));
fprintf(fp,"--------------------------------\n");
fprintf(fp, "        N ATOM  RES RESN CH. CHARGE  COUL_NRG  SOLV_NRG  BORN_NRG   POL_NRG  SASA_NRG   TOTAL_NRG (kJ/mol)\n"); 
for(i=0; i<sys.n_atoms; i++)
   fprintf(fp, "NRG %5i %4s %4s %4i %c  %7.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf\n", 
               sys.atoms[i].at_n,
               sys.atoms[i].at_name,
               sys.atoms[i].res_name,
               sys.atoms[i].res_n,
               sys.atoms[i].chain,
               sys.atoms[i].charge,
               nrg.coul_a[i],
               nrg.solv_a[i],
               nrg.born_a[i],
               nrg.pol_a[i],
               nrg.area_a[i],
               nrg.tot_a[i]
               );
fclose(fp);
        printf("--------------------------------\n");
}

void out_grid(struct System sys, struct Atom_grid atom_grid, double *gbr6, struct Neighbour *neighbours, struct Const constants, struct Flag_par flag_par)
{
	double xmin,xmax,ymin,ymax,zmin,zmax,t1;
	double hx, hy, hz;
	double f, d, *v;
        int i, j, k, l, ii, jj, kk, m;
	int nx, ny, nz,il,jl,kl;
	double ***out;
	bool ***isin;
	char buf[1024];
	FILE *fp;
	nx = flag_par.nx;
	ny = flag_par.ny;
	nz = flag_par.nz;
v = (double *) calloc((size_t) 3, sizeof(double));
if(flag_par.grid)
{
if(nx > 0 && ny > 0 && nz > 0)
{
if(1.0/constants.kd > 10.0) t1 = 0.1; else t1 = constants.kd;
xmin = atom_grid.x_min - 3.0/t1;
xmax = atom_grid.x_max + 3.0/t1;
ymin = atom_grid.y_min - 3.0/t1;
ymax = atom_grid.y_max + 3.0/t1;
zmin = atom_grid.z_min - 3.0/t1;
zmax = atom_grid.z_max + 3.0/t1;
hx = (xmax - xmin)/(double) nx;
hy = (ymax - ymin)/(double) ny;
hz = (zmax - zmin)/(double) nz;
}
else
{
xmin = (atom_grid.x_min + atom_grid.x_max)/2.0 - (double) nx * flag_par.mesh/2.0;
ymin = (atom_grid.y_min + atom_grid.y_max)/2.0 - (double) ny * flag_par.mesh/2.0;
zmin = (atom_grid.z_min + atom_grid.z_max)/2.0 - (double) nz * flag_par.mesh/2.0;
xmin = (atom_grid.x_min + atom_grid.x_max)/2.0 + (double) nx * flag_par.mesh/2.0;
ymin = (atom_grid.y_min + atom_grid.y_max)/2.0 + (double) ny * flag_par.mesh/2.0;
zmin = (atom_grid.z_min + atom_grid.z_max)/2.0 - (double) nz * flag_par.mesh/2.0;
hx = hy = hz = flag_par.mesh;
}

out = (double ***) calloc((size_t) nx, sizeof(double **));
for(i = 0; i < nx; i++)
{
out[i] = (double **) calloc((size_t) ny, sizeof(double *));
for(j = 0; j < ny; j++)
out[i][j] = (double *) calloc((size_t) nz, sizeof(double));
}
isin = (bool ***) calloc((size_t) nx, sizeof(bool **));
for(i = 0; i < nx; i++)
{
isin[i] = (bool **) calloc((size_t) ny, sizeof(bool *));
for(j = 0; j < ny; j++)
isin[i][j] = (bool *) calloc((size_t) nz, sizeof(bool));
}

strcpy(buf,flag_par.file_out);
strcat(buf,".dx");
fp = fopen(buf,"w");
printf(".....I am writing the potential in DX format in\n%s\n",buf);
fprintf(fp,"# Comments\n\
object 1 class gridpositions counts %i %i %i\n\
origin %f %f %f\n\
delta %f 0.0 0.0\n\
delta 0.0 %f 0.0\n\
delta 0.0 0.0 %f\n\
object 2 class gridconnections counts %i %i %i\n\
object 3 class array type double rank 0 items %i data follows\n",
nx, ny, nz, xmin, ymin, zmin,hx,hy,hz,nx,ny,nz, nx*ny*nz );
if( (flag_par.cutoff > 3.0/t1 - hx || flag_par.cutoff > 3.0/t1 - hy || flag_par.cutoff > 3.0/t1 - hz) && (nx == -1 || ny == -1 || nz == -1)) 
{
printf("Reset cutoff for grid calculation to %lf...\n", 3.0/t1);
il = (int) (3.0/t1)/hx - 1;
jl = (int) (3.0/t1)/hy - 1;
kl = (int) (3.0/t1)/hz - 1;
}
else
{
il = (int) (flag_par.cutoff)/hx + 1;
jl = (int) (flag_par.cutoff)/hy + 1;
kl = (int) (flag_par.cutoff)/hz + 1;
}
for(m=0; m<sys.n_atoms; m++)
{
i = (int) ((sys.atoms[m].coor[0] - xmin) / hx);
j = (int) ((sys.atoms[m].coor[1] - ymin) / hy);
k = (int) ((sys.atoms[m].coor[2] - zmin) / hz);
for(ii = -il; ii <= il ; ii++)
for(jj = -jl; jj <= jl ; jj++)
for(kk = -kl; kk <= kl ; kk++)
if(!isin[i+ii][j+jj][k+kk])
{
v[0] = xmin + (double) (i+ii) * hx;
v[1] = ymin + (double) (j+jj) * hy;
v[2] = zmin + (double) (k+kk) * hz;
d = distv( sys.atoms[m].coor, v);
f = fgb(gbr6[i],flag_par.test_charge_radius,  flag_par.probe_radius,  flag_par.salt_radius,  d, flag_par.pdie,  flag_par.sdie,  constants.kd, flag_par.kp);
out[i+ii][j+jj][k+kk] = out[i+ii][j+jj][k+kk] + sys.atoms[m].charge/(d *flag_par.pdie); 
out[i+ii][j+jj][k+kk] = out[i+ii][j+jj][k+kk] + sys.atoms[m].charge*f; 
if (d < (sys.atoms[m].radius + flag_par.probe_radius)) 
{
out[i+ii][j+jj][k+kk]=0.0;
isin[i+ii][j+jj][k+kk]=1;
}
}
}

for(i = 0,l=1; i < nx; i++)
for(j = 0; j < ny; j++)
for(k = 0; k < nz; k++)
{
out[i][j][k] = constants.k_el * out[i][j][k] /1000;
fprintf(fp, "%f ", out[i][j][k]);
if(l%3==0) fprintf(fp,"\n");
l++;
}
if( (l+1)%3 != 2) fprintf(fp,"\n");
fprintf(fp,"attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n");
fclose(fp);
        printf("--------------------------------\n");
}

}


void get_grid_parameters(struct System sistema, double probe_radius, 
		double cutoff, struct Atom_grid *atom_grid, char type)
{

	int i;
	double safe = 0.1, max_vdw_ctc = 0.0, r_min, min_vol_atom;

        for(i=0; i< sistema.n_atoms; i++)
        if((2.0 * sistema.atoms[i].radius) > max_vdw_ctc)
            max_vdw_ctc = 2.0 * sistema.atoms[i].radius;
        r_min = max_vdw_ctc/2.0; 
        for(i=0; i< sistema.n_atoms; i++)
         if((sistema.atoms[i].radius) < r_min)
          r_min = sistema.atoms[i].radius;
        min_vol_atom = 4 * M_PI * r_min * r_min;

if(type='v')
{
if (probe_radius >= (cutoff/2.0) ) 
	(*atom_grid).mesh = max_vdw_ctc + 2*probe_radius + safe;
else
	(*atom_grid).mesh = max_vdw_ctc + cutoff + safe;
}

if(type=='d')
(*atom_grid).mesh = cutoff + safe;

for (i=0; i< sistema.n_atoms; i++)
{
	 if(i == 0)
	     {
	       (*atom_grid).x_min = (*atom_grid).x_max = sistema.atoms[i].coor[0];
	       (*atom_grid).y_min = (*atom_grid).y_max = sistema.atoms[i].coor[1];
               (*atom_grid).z_min = (*atom_grid).z_max = sistema.atoms[i].coor[2];
              }
               if(sistema.atoms[i].coor[0] > (*atom_grid).x_max) (*atom_grid).x_max = sistema.atoms[i].coor[0];
               if(sistema.atoms[i].coor[1] > (*atom_grid).y_max) (*atom_grid).y_max = sistema.atoms[i].coor[1];
               if(sistema.atoms[i].coor[2] > (*atom_grid).z_max) (*atom_grid).z_max = sistema.atoms[i].coor[2];
               if(sistema.atoms[i].coor[0] < (*atom_grid).x_min) (*atom_grid).x_min = sistema.atoms[i].coor[0];
               if(sistema.atoms[i].coor[1] < (*atom_grid).y_min) (*atom_grid).y_min = sistema.atoms[i].coor[1];
               if(sistema.atoms[i].coor[2] < (*atom_grid).z_min) (*atom_grid).z_min = sistema.atoms[i].coor[2];

}
(*atom_grid).x_min = (*atom_grid).x_min - safe;
(*atom_grid).y_min = (*atom_grid).y_min - safe;
(*atom_grid).z_min = (*atom_grid).z_min - safe;
(*atom_grid).x_max = (*atom_grid).x_max + safe;
(*atom_grid).y_max = (*atom_grid).y_max + safe;
(*atom_grid).z_max = (*atom_grid).z_max + safe;

(*atom_grid).max_atom_node = (int) ((*atom_grid).mesh * (*atom_grid).mesh * 
		                  (*atom_grid).mesh / (double) min_vol_atom ) + 1;


            (*atom_grid).grid_X = 
	    (int) (((*atom_grid).x_max - (*atom_grid).x_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_Y = 
	    (int) (((*atom_grid).y_max - (*atom_grid).y_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_Z = 
	    (int) (((*atom_grid).z_max - (*atom_grid).z_min) / (*atom_grid).mesh ) + 1;
            
	    (*atom_grid).grid_size = 
	    (*atom_grid).grid_X * (*atom_grid).grid_Y * (*atom_grid).grid_Z;


if((*atom_grid).PBC_ON == 1)
{
(*atom_grid).dx = (*atom_grid).pbcx / (double) (*atom_grid).grid_X;
(*atom_grid).dy = (*atom_grid).pbcy / (double) (*atom_grid).grid_Y;
(*atom_grid).dz = (*atom_grid).pbcz / (double) (*atom_grid).grid_Z;
if( 
((*atom_grid).x_max - (*atom_grid).x_min) > (*atom_grid).pbcx ||
((*atom_grid).y_max - (*atom_grid).y_min) > (*atom_grid).pbcy ||
((*atom_grid).z_max - (*atom_grid).z_min) > (*atom_grid).pbcz 
)
{
printf("Warning: PBC less than system dimension...\ncontinuing....\n");
}
}
}

void atoms_on_grid(struct Atom_grid *atom_grid, struct System sistema)
{
	int i,j,k,l,n, index;
        double dx,dy,dz, safe = 0.1;




for(i=0; i < (*atom_grid).grid_size; i++)
(*atom_grid).n_atom_node[i]=0;
if((*atom_grid).PBC_ON != 1)
for(n=0; n< sistema.n_atoms; n++)
{
	i = (int) ((sistema.atoms[n].coor[0] - (*atom_grid).x_min)/((*atom_grid).mesh ));
	j = (int) ((sistema.atoms[n].coor[1] - (*atom_grid).y_min)/((*atom_grid).mesh ));
	k = (int) ((sistema.atoms[n].coor[2] - (*atom_grid).z_min)/((*atom_grid).mesh ));
if(i >=  (*atom_grid).grid_X || j >=  (*atom_grid).grid_Y || k >= (*atom_grid).grid_Z || i < 0 || j < 0 || k < 0)
{
printf("%i %i %i vs %i %i %i\n", i, j, k, (*atom_grid).grid_X,(*atom_grid).grid_Y,(*atom_grid).grid_Z);
printf("%lf %lf %lf - %lf %lf %lf - %lf %lf %lf\n", 
	(*atom_grid).x_min, sistema.atoms[n].coor[0], (*atom_grid).x_max,
	(*atom_grid).y_min, sistema.atoms[n].coor[1], (*atom_grid).y_max,
	(*atom_grid).z_min, sistema.atoms[n].coor[2], (*atom_grid).z_max
	);
exit(0);
}
	index = i * (*atom_grid).grid_Y *(*atom_grid).grid_Z + j * (*atom_grid).grid_Z + k;
        (*atom_grid).atom_node[index][(*atom_grid).n_atom_node[index]] = n;
        (*atom_grid).n_atom_node[index] = (*atom_grid).n_atom_node[index] + 1;
if( (*atom_grid).n_atom_node[index] >= (*atom_grid).max_atom_node) 
{
(*atom_grid).max_atom_node =  (int) ((double) (*atom_grid).max_atom_node * sqrt(2));
for(l = 0; l< (*atom_grid).grid_size; l++)
(*atom_grid).atom_node[l] = realloc((*atom_grid).atom_node[l], (size_t) (*atom_grid).max_atom_node * sizeof(int));
}


}
else
for(n=0; n< sistema.n_atoms; n++)
{
	i = ((int) ((sistema.atoms[n].coor[0] - (*atom_grid).x_min)/(*atom_grid).dx))%(*atom_grid).grid_X ;
	j = ((int) ((sistema.atoms[n].coor[1] - (*atom_grid).y_min)/(*atom_grid).dy))%(*atom_grid).grid_Y ;
	k = ((int) ((sistema.atoms[n].coor[2] - (*atom_grid).z_min)/(*atom_grid).dz))%(*atom_grid).grid_Z ;
	index = i * (*atom_grid).grid_Y *(*atom_grid).grid_Z + j * (*atom_grid).grid_Z + k;

	(*atom_grid).atom_node[index][(*atom_grid).n_atom_node[index]] = n;
	(*atom_grid).n_atom_node[index]++; 
}
j=0;





for(l = 0; l< (*atom_grid).grid_size; l++)
(*atom_grid).atom_node[l] = realloc((*atom_grid).atom_node[l], (*atom_grid).n_atom_node[l] * sizeof(int));
}


void grid_neighbour(struct Atom_grid atom_grid, struct System sistema, 
double cutoff, struct Neighbour *neighbours, char type)
{
int i,j,k,r,s,t,index1,index2,m,n;
int i1,i2;
double x1,y1,z1,x2,y2,z2,d,tt;
double DX, DY, DZ;
for (i=0; i<sistema.n_atoms; i++) neighbours[i].n_neighbours = 0;


for(i = 0; i< atom_grid.grid_X; i++)
for(j = 0; j< atom_grid.grid_Y; j++)
for(k = 0; k< atom_grid.grid_Z; k++)
{
index1 = i*atom_grid.grid_Y *atom_grid.grid_Z + j * atom_grid.grid_Z + k;
for(r = -1; r<= 1; r++)
for(s = -1; s<= 1; s++)
for(t = -1; t<= 1; t++)
if(atom_grid.PBC_ON != 1)
{
index2 = (i+r)*atom_grid.grid_Y *atom_grid.grid_Z + 
	 (j+s)*atom_grid.grid_Z + k+t;
	if( ((i+r) < atom_grid.grid_X ) && ((j+s) < atom_grid.grid_Y) 
		&& ((k+t) < atom_grid.grid_Z) && 
	     ((i+r) >= 0 ) && ((j+s) >= 0) && ((k+t) >= 0) )
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
	{
		i1 = atom_grid.atom_node[index1][n];
		x1 = sistema.atoms[i1].coor[0]; 
		y1 = sistema.atoms[i1].coor[1]; 
		z1 = sistema.atoms[i1].coor[2]; 
		for(m=0; m< atom_grid.n_atom_node[index2]; m++)
		{
		i2 = atom_grid.atom_node[index2][m];
		x2 = sistema.atoms[i2].coor[0]; 
		y2 = sistema.atoms[i2].coor[1]; 
		z2 = sistema.atoms[i2].coor[2];
        if(type = 'v')
	d =  (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff)*
             (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff);
        else if(type = 'd') d = cutoff*cutoff;

	if( ((x1-x2)*(x1-x2)) <= d) 
	if( ((y1-y2)*(y1-y2)) <= d) 
	if( ((z1-z2)*(z1-z2)) <= d)
        {
	tt = ((x1-x2)*(x1-x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2)); 
        if(tt <= d)
        if(i1 > i2)
        {
        neighbours[i1].list[neighbours[i1].n_neighbours] = i2;
        neighbours[i2].list[neighbours[i2].n_neighbours] = i1;
        neighbours[i1].d[neighbours[i1].n_neighbours] = sqrt(tt);
        neighbours[i2].d[neighbours[i2].n_neighbours] = sqrt(tt);
        neighbours[i1].n_neighbours++;
        neighbours[i2].n_neighbours++;
        if(neighbours[i1].n_neighbours >= neighbours[i1].max_n)
         {
         neighbours[i1].max_n = neighbours[i1].max_n * 2;
         neighbours[i1].list = realloc(neighbours[i1].list, (size_t) neighbours[i1].max_n * sizeof(int));
         neighbours[i1].d = realloc(neighbours[i1].d, neighbours[i1].max_n * sizeof(double));
         }
        if(neighbours[i2].n_neighbours >= neighbours[i2].max_n)
         {
         neighbours[i2].max_n = neighbours[i2].max_n * 2;
         neighbours[i2].list = realloc(neighbours[i2].list, (size_t) neighbours[i2].max_n * sizeof(int));
         neighbours[i2].d = realloc(neighbours[i2].d, neighbours[i2].max_n * sizeof(double));
         }
	}
        }
        }
}
else if(atom_grid.PBC_ON == 1) 
{
DX = atom_grid.pbcx;
DY = atom_grid.pbcy;
DZ = atom_grid.pbcz;
index2 = ((i+r+atom_grid.grid_X)%atom_grid.grid_X)*atom_grid.grid_Y *atom_grid.grid_Z + 
	 ((j+s+atom_grid.grid_Y)%atom_grid.grid_Y)*atom_grid.grid_Z + (k+t+atom_grid.grid_Z)%atom_grid.grid_Z;
if((atom_grid.grid_X > 2 || ((atom_grid.grid_X == 2) && r >= 0) || ((atom_grid.grid_X == 1) && r == 0)) 
&& (atom_grid.grid_Y > 2 || ((atom_grid.grid_Y == 2) && s >= 0) || ((atom_grid.grid_Y == 1) && s == 0))
&& (atom_grid.grid_Z > 2 || ((atom_grid.grid_Z == 2) && t >= 0) || ((atom_grid.grid_Z == 1) && t == 0))
)
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
	{
		i1 = atom_grid.atom_node[index1][n];
		x1 = sistema.atoms[i1].coor[0]; 
		y1 = sistema.atoms[i1].coor[1]; 
		z1 = sistema.atoms[i1].coor[2]; 
		for(m=0; m< atom_grid.n_atom_node[index2]; m++)
		{
		i2 = atom_grid.atom_node[index2][m];
		x2 = sistema.atoms[i2].coor[0]; 
		y2 = sistema.atoms[i2].coor[1]; 
		z2 = sistema.atoms[i2].coor[2];
        if(type = 'v')
	d =  (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff)*
             (sistema.atoms[i1].radius + sistema.atoms[i2].radius + cutoff);
        else if(type = 'd') d = cutoff*cutoff;

        x2 = (fabs(x1-x2) < fabs(x1-x2+DX)) ? (x2) : (x2-DX);
        x2 = (fabs(x1-x2) < fabs(x1-x2-DX)) ? (x2) : (x2+DX);
        y2 = (fabs(y1-y2) < fabs(y1-y2+DY)) ? (y2) : (y2-DY);
        y2 = (fabs(y1-y2) < fabs(y1-y2-DY)) ? (y2) : (y2+DY);
        z2 = (fabs(z1-z2) < fabs(z1-z2+DZ)) ? (z2) : (z2-DZ);
        z2 = (fabs(z1-z2) < fabs(z1-z2-DZ)) ? (z2) : (z2+DZ);
 
	if( ((x1-x2)*(x1-x2)) <= d) 
	if( ((y1-y2)*(y1-y2)) <= d) 
	if( ((z1-z2)*(z1-z2)) <= d)
        {
	tt = ((x1-x2)*(x1-x2)+(y1 - y2)*(y1 - y2)+(z1 - z2)*(z1 - z2));
        if(tt <= d)
        if(i1 > i2)
        {
        neighbours[i1].list[neighbours[i1].n_neighbours] = i2;
        neighbours[i2].list[neighbours[i2].n_neighbours] = i1;
        neighbours[i1].d[neighbours[i1].n_neighbours] = sqrt(tt);
        neighbours[i2].d[neighbours[i2].n_neighbours] = sqrt(tt);
        neighbours[i1].n_neighbours++;
        neighbours[i2].n_neighbours++;

        if(neighbours[i1].n_neighbours >= neighbours[i1].max_n)
         {
         neighbours[i1].max_n = neighbours[i1].max_n * 2;
         neighbours[i1].list = realloc(neighbours[i1].list, (size_t) neighbours[i1].max_n * sizeof(int));
         neighbours[i1].d = realloc(neighbours[i1].d, neighbours[i1].max_n * sizeof(double));
         }
        if(neighbours[i2].n_neighbours >= neighbours[i2].max_n)
         {
         neighbours[i2].max_n = neighbours[i2].max_n * 2;
         neighbours[i2].list = realloc(neighbours[i2].list, (size_t) neighbours[i2].max_n * sizeof(int));
         neighbours[i2].d = realloc(neighbours[i2].d, neighbours[i2].max_n * sizeof(double));
         }
        }
        }
	}
                 }
        }
}

}
}

void system_area_alloc(struct System system, struct System_area *system_area)
{
(*system_area).atoms=(struct Area *) calloc((size_t) system.n_atoms, sizeof(struct Area));
(*system_area).residues=(struct Area *) calloc((size_t) system.n_residues, sizeof(struct Area));
(*system_area).chains=(struct Area *) calloc((size_t) system.n_chains, sizeof(struct Area));
(*system_area).segments=(struct Area *) calloc((size_t) system.n_segments, sizeof(struct Area));
(*system_area).models=(struct Area *) calloc((size_t) system.n_models, sizeof(struct Area));
}

void system_area_free(struct System_area *system_area)
{
free((*system_area).atoms);
free((*system_area).residues);
free((*system_area).chains);
free((*system_area).segments);
free((*system_area).models);
}


void grid_srf_dens(
                struct System system,
                struct Atom_grid atom_grid,
                struct Flag_par flag_par,
                struct Srf *srf)
{
int i,j,k,kk,r,s,t,m,n,i1,i2,ii,nsrfpt, ir, iv, n_radii, i_srf=0;
int index1, index2, *n_srf_pt;
double acc, acc_self_chain, acc_self_res;
double phi, psi, area_sas, da = 0.0, area_ctc, d, tmpf;
double r_min, r_max, dr;
double **srfxr, **srfyr, **srfzr,*rad, dx, dist2, dist;
double *srfx, *srfy, *srfz;
double *sas_r, *sas_c;
double x1, y1, z1, x2, y2, z2, x, y, z;
double fdt,ft,fds,fs,fu;
r_max = 0.0;
r_min = DBL_MAX;
        for(i=0; i< system.n_atoms; i++)
        {
        if((system.atoms[i].radius) > r_max)
            r_max = system.atoms[i].radius ;
        if((system.atoms[i].radius ) < r_min)
            r_min = system.atoms[i].radius ;

        }
dr = 0.01;
n_radii = (int) floor(0.000001 + (r_max - r_min)/dr) + 1;
        srfxr = (double **) calloc((size_t) n_radii, sizeof(double *));
        srfyr = (double **) calloc((size_t) n_radii, sizeof(double *));
        srfzr = (double **) calloc((size_t) n_radii, sizeof(double *));
        sas_r = (double *) calloc((size_t) n_radii, sizeof(double));
        sas_c = (double *) calloc((size_t) n_radii, sizeof(double));
        rad   = (double *) calloc((size_t) n_radii, sizeof(double));
        n_srf_pt  = (int *) calloc((size_t) n_radii, sizeof(int));
for(i = 0; i < system.n_atoms; i++) 
(*srf).at_srf[i].n_isrf  = 0;
for(i = 0; i < n_radii; i++) 
{
	n_srf_pt[i] = 4.0 * M_PI * (0.0001 + r_min + dr * (double) i ) * (0.0001 + r_min + dr * (double) i )/flag_par.msms_area;
        srfxr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
        srfyr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
        srfzr[i] = (double *) calloc(n_srf_pt[i], sizeof(double));
}

for(x = r_min; x<= r_max ; x = x + dr)
{
i = (int) floor(0.0001 + (x - r_min)/dr);

rad[i] = x; 
fdt = M_PI * (3-sqrt(5));
ft=0.0;
fds = 2.0/(double) n_srf_pt[i];
fs = 1 - fds/2.0;
for(j = 0; j< n_srf_pt[i]; j++)
{
fu = sqrt(1-fs*fs);
srfxr[i][j] =  (x + flag_par.probe_radius ) * cos(ft)*fu;
srfyr[i][j] =  (x + flag_par.probe_radius ) * sin(ft)*fu;
srfzr[i][j] =  (x + flag_par.probe_radius ) * fs;
fs = fs -fds;
ft = ft + fdt;
}
sas_r[i] = (x + flag_par.probe_radius) * (x + flag_par.probe_radius) *4*M_PI/(double) n_srf_pt[i];
sas_c[i] = (x) * (x) *4*M_PI/(double) n_srf_pt[i];
}
for(i = 0; i< atom_grid.grid_X; i++)
for(j = 0; j< atom_grid.grid_Y; j++)
for(k = 0; k< atom_grid.grid_Z; k++)
{
     index1 = i*atom_grid.grid_Y *atom_grid.grid_Z + j * atom_grid.grid_Z + k;
        for(n=0; n< atom_grid.n_atom_node[index1]; n++)
          {
                i1 = atom_grid.atom_node[index1][n];
                if ((system.atoms[i1].radius >= r_min) && (system.atoms[i1].radius <= r_max)) 
                {
                ir = (int) floor( (system.atoms[i1].radius - r_min)/dr + 0.0001);
                area_sas = sas_r[ir];
                srfx = srfxr[ir]; srfy = srfyr[ir]; srfz = srfzr[ir];
                }
                else printf("I don't know ATOM: element %s: %s in res %s \n", system.atoms[i1].element,  system.atoms[i1].at_name, system.atoms[i1].res_name);

        if(system.atoms[i1].radius < r_min || system.atoms[i1].radius > r_max)
            {
            printf("%s %s %i %c has radius %11.3f... reset to 1.5\n",
            system.atoms[i1].at_name,
            system.atoms[i1].res_name,
            system.atoms[i1].res_n,
            system.atoms[i1].chain,
            system.atoms[i1].radius);
            system.atoms[i1].radius = 1.5;
            }
        for(ii = 0; ii<n_srf_pt[ir]; ii++)
        {
        x1 = system.atoms[i1].coor[0] + srfx[ii];
        y1 = system.atoms[i1].coor[1] + srfy[ii];
        z1 = system.atoms[i1].coor[2] + srfz[ii];

        acc_self_res = acc_self_chain = acc = 1.0;
        for(r = -1; r<= 1; r++)
        for(s = -1; s<= 1; s++)
        for(t = -1; t<= 1; t++)
        {
        index2 = (i+r)*atom_grid.grid_Y * atom_grid.grid_Z + (j+s)*atom_grid.grid_Z + k+t;
        if(
           ((i+r) < atom_grid.grid_X ) && ((j+s) < atom_grid.grid_Y) &&
           ((k+t) < atom_grid.grid_Z) && ((i+r) >= 0 ) && ((j+s) >= 0) &&
           ((k+t) >= 0)
          )
        for(m=0; m< atom_grid.n_atom_node[index2]; m++)
        if(atom_grid.atom_node[index2][m] != i1)
        {
        i2 = atom_grid.atom_node[index2][m];

        x2 = system.atoms[i2].coor[0];
        y2 = system.atoms[i2].coor[1];
        z2 = system.atoms[i2].coor[2];
        d = system.atoms[i2].radius+flag_par.probe_radius;
        x = x1-x2;
        if(fabs(x) < d)
         {
          y = y1-y2;
           if(fabs(y) < d)
            {
             z = z1-z2;
              if(fabs(z) < d)
               if(i1!=i2)
               {
               dist2 = x*x + y*y + z*z;
               if(dist2 < (d*d))
                  {
                  acc = 0;

                  if(
                     ((system.atoms[i2].chain == system.atoms[i1].chain) && 
                      !(strcmp(system.atoms[i1].segid,system.atoms[i2].segid)))
                    )
                    {
                    acc_self_chain = 0; 
                    if(system.atoms[i2].res_n == system.atoms[i1].res_n)
                    if(system.atoms[i2].res_ins == system.atoms[i1].res_ins)
                    acc_self_res = 0; 
                    }
                  }
              }
            }
          }
        if(acc == 0) 
        m=atom_grid.n_atom_node[index2];
         }
	}
        if(acc > 0.0)
        {
        if(i_srf >= ((*srf).alloc_srf_pt)) {
printf("Attempting to write more surface points than allocated, I will reallocate surface points...\n");
         (*srf).alloc_srf_pt = (int) ((double) (*srf).alloc_srf_pt * sqrt(2));
         (*srf).srf_pt = realloc((*srf).srf_pt, sizeof(struct Srf_pt) * (*srf).alloc_srf_pt);
         for(kk=i_srf; kk< (*srf).alloc_srf_pt ; kk++)
            {
            (*srf).srf_pt[kk].r = (double *) calloc((size_t) 3, sizeof(double));
            (*srf).srf_pt[kk].vec = (double *) calloc((size_t) 3, sizeof(double));
            }
         if((*srf).srf_pt == NULL) printf("Could not allocate %i srf points\n...exiting\n", (*srf).alloc_srf_pt ); exit(0);
         }
        (*srf).srf_pt[i_srf].a = area_sas  * acc ;
        (*srf).srf_pt[i_srf].r[0] = x1;
        (*srf).srf_pt[i_srf].r[1] = y1;
        (*srf).srf_pt[i_srf].r[2] = z1;
        (*srf).srf_pt[i_srf].vec[0] = srfx[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].vec[1] = srfy[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].vec[2] = srfz[ii]/(rad[ir] + flag_par.probe_radius);
        (*srf).srf_pt[i_srf].at = i1;
        (*srf).at_srf[i1].isrf[ (*srf).at_srf[i1].n_isrf++] = i_srf;
        i_srf++;
        }
   }
  }
}
for(kk = i_srf; kk<(*srf).alloc_srf_pt; kk++)
{
free((*srf).srf_pt[kk].r);
free((*srf).srf_pt[kk].vec);
}
(*srf).n_srf_pt_tot = i_srf;
(*srf).alloc_srf_pt = i_srf;
         (*srf).srf_pt = (struct Srf_pt *) realloc((*srf).srf_pt, (size_t) (sizeof(struct Srf_pt) * (*srf).alloc_srf_pt));

for(i=0,(*srf).area=0; i<(*srf).n_srf_pt_tot; i++)
(*srf).area=(*srf).area + (*srf).srf_pt[i].a;

for(i = 0; i < n_radii; i++) 
{
        free(srfxr[i]); 
        free(srfyr[i]); 
        free(srfzr[i]); 
}
        free(srfxr); 
        free(srfyr); 
        free(srfzr); 
        free(sas_r); 


}

void print_grid_info(struct Atom_grid atom_grid)
{
printf("grid (X,Y,Z): %i %i %i -- grid size: %i -- grid mesh: %lf\n", atom_grid.grid_X, atom_grid.grid_Y, atom_grid.grid_Z, atom_grid.grid_size, atom_grid.mesh);
if(atom_grid.PBC_ON)
printf("x: %lf -- %lf    dx: %lf\ny: %lf -- %lf    dy: %lf\nz: %lf -- %lf    dz: %lf\n",
atom_grid.x_min, atom_grid.x_max, atom_grid.dx,
atom_grid.y_min, atom_grid.y_max, atom_grid.dy,
atom_grid.z_min, atom_grid.z_max, atom_grid.dz);
else 
printf("x: %lf -- %lf\ny: %lf -- %lf\nz: %lf -- %lf\n",
atom_grid.x_min, atom_grid.x_max,
atom_grid.y_min, atom_grid.y_max,
atom_grid.z_min, atom_grid.z_max);
}

void srf_free(struct Srf *srf, int n_atoms)
{
int i;
for(i=0;i<(*srf).alloc_srf_pt;i++)
{
free((*srf).srf_pt[i].r);
free((*srf).srf_pt[i].vec);
} 
free((*srf).srf_pt);
for(i=0;i<n_atoms;i++)
free((*srf).at_srf[i].isrf);
free((*srf).at_srf);
}

int is_within(struct Atom *atoms,int n_atoms,struct Atom probe_atom, double cutoff)
{
static int i,j;
for(i=0,j=0; j< n_atoms; j++)
if(distv(atoms[j].coor, probe_atom.coor) <= cutoff)
{
i=1;
break;
}
return i;
}

void print_summary(struct Summary sum, struct Flag_par flag_par)
{
char buf[2048];
FILE *fpsum;

sprintf(buf,"%s.summary", flag_par.file_out);

fpsum = file_open(buf,"w");
        fprintf(fpsum,"--------------------------------\nRUN PARAMETERS:\n--------------------------------\n");
        fprintf(fpsum,"pdb file: %s\n", flag_par.file_in_pdb);
        fprintf(fpsum,"output file: %s\n", flag_par.file_out);
        fprintf(fpsum,"probe_radius: %lf\n", flag_par.probe_radius);
        fprintf(fpsum,"minimum atomic radius: %lf\n", flag_par.min_radius);
        fprintf(fpsum,"ionic strength (M): %lf\n", flag_par.ions);
        fprintf(fpsum,"salt_radius (A): %lf\n", flag_par.salt_radius);
        fprintf(fpsum,"inner dielectric constant: %lf\n", flag_par.pdie);
        fprintf(fpsum,"outer dielectric constant: %lf\n", flag_par.sdie);
        fprintf(fpsum,"temperature (K): %lf\n", flag_par.temp);
        fprintf(fpsum,"cutoff for GBR6 calculation (A): %lf\n", flag_par.cutoff);
        fprintf(fpsum,"cutoff for electrostatic patches (A): %lf\n", flag_par.cutoff2);
        fprintf(fpsum,"cutoff for contact definition between vdw surfs. (A): %lf\n", flag_par.cutoff_short);
        fprintf(fpsum,"gamma (kJ/A^2): %lf\n", flag_par.gamma);
        if(flag_par.srfpot || flag_par.grid)
        {
        fprintf(fpsum,"test charge radius: %lf\n", flag_par.test_charge_radius);
        }
        if(flag_par.msms )
        {
        fprintf(fpsum,"base file name for MSMS: %s\n", flag_par.file_msms);
        }
        if(flag_par.ns )
        {
        fprintf(fpsum,"base filename for NanoShaper: %s\n", flag_par.file_ns);
        }
        fprintf(fpsum,"area (A^2) per point on surface: %lf\n", flag_par.msms_area);
        fprintf(fpsum,"--------------------------------\n\n");
        fprintf(fpsum,"Charge at interface\n");
        fprintf(fpsum,"--------------------------------\n");
fprintf(fpsum,"Total charge on complex interface region: %11.3f e\n", sum.charge);
fprintf(fpsum,"Total charge on molecule 1 interface region: %11.3f e\n", sum.charge_mol_1);
fprintf(fpsum,"Total charge on molecule 2 interface region: %11.3f e\n", sum.charge_mol_2);
fprintf(fpsum,"Sum of products of contacting atoms charges at interface: %11.3f  e^2\n", sum.ccsum);
fprintf(fpsum,"Pearson correlation of contacting atoms charges at interface: %6.3f\n", sum.ccpc);
        fprintf(fpsum,"--------------------------------\n\n");

fprintf(fpsum,"Interaction energy\n");
fprintf(fpsum,"--------------------------------\n");
fprintf(fpsum,"Total Coulomb   energy: %14.4lf (kJ/mol)\n", sum.coul_nrg);
fprintf(fpsum,"Total solvation energy: %14.4lf (kJ/mol)\n", sum.solv_nrg);
fprintf(fpsum," --Born self    energy: %14.4lf (kJ/mol)\n", sum.born_nrg);
fprintf(fpsum," --Coulomb solv energy: %14.4lf (kJ/mol)\n", sum.coul_solv_nrg);
fprintf(fpsum," --sasa         energy: %14.4lf (kJ/mol)\n", sum.sasa_nrg);
fprintf(fpsum,"Total           energy: %14.4lf (kJ/mol)\n", (sum.coul_nrg + sum.solv_nrg));
fprintf(fpsum,"--------------------------------\n\n");

        fprintf(fpsum,"Surface buried at interface\n");
        fprintf(fpsum,"--------------------------------\n");
fprintf(fpsum,"Solvent accessible surface area buried at interface: %11.3f A^2\n", sum.sas);
if(flag_par.msms || flag_par.ns)
fprintf(fpsum,"Solvent excluded surface area buried at interface: %11.3f A^2\n", sum.ses);
fprintf(fpsum,"--------------------------------\n\n");
if(flag_par.msms || flag_par.ns)
{
        fprintf(fpsum,"Electrostatic potential at interface\n");
        fprintf(fpsum,"--------------------------------\n");
fprintf(fpsum,"Sum of products of contacting surface points potentials at interface: %11.3f (kJ/ (e mol))^2\n", sum.ppsum);
fprintf(fpsum,"Pearson correlation of contacting surface points potentials at interface: %6.3f\n", sum.pppc);
        fprintf(fpsum,"--------------------------------\n\n");
}
else
{
        fprintf(fpsum,"No information on Electrostatic potential at interface\n");
        fprintf(fpsum,"because nor NanoShaper or msms was used to generate the solvent excluded surface\n");
}

        printf("--------------------------------\nRUN PARAMETERS:\n--------------------------------\n");
        printf("pdb file: %s\n", flag_par.file_in_pdb);
        printf("output file: %s\n", flag_par.file_out);
        printf("probe_radius: %lf\n", flag_par.probe_radius);
        printf("minimum atomic radius: %lf\n", flag_par.min_radius);
        printf("ionic strength (M): %lf\n", flag_par.ions);
        printf("salt_radius (A): %lf\n", flag_par.salt_radius);
        printf("inner dielectric constant: %lf\n", flag_par.pdie);
        printf("outer dielectric constant: %lf\n", flag_par.sdie);
        printf("temperature (K): %lf\n", flag_par.temp);
        printf("cutoff for GBR6 calculation (A): %lf\n", flag_par.cutoff);
        printf("cutoff for electrostatic patches (A): %lf\n", flag_par.cutoff2);
        printf("cutoff for contact definition between vdw surfs. (A): %lf\n", flag_par.cutoff_short);
        printf("gamma (kJ/A^2): %lf\n", flag_par.gamma);
        if(flag_par.srfpot || flag_par.grid)
        {
        printf("test charge radius: %lf\n", flag_par.test_charge_radius);
        }
        if(flag_par.msms )
        {
        printf("base file name for MSMS: %s\n", flag_par.file_msms);
        }
        if(flag_par.ns )
        {
        printf("base filename for NanoShaper: %s\n", flag_par.file_ns);
        }
        printf("area (A^2) per point on surface: %lf\n", flag_par.msms_area);
        printf("--------------------------------\n\n");
        printf("Charge at interface\n");
        printf("--------------------------------\n");
printf("Total charge on complex interface region: %11.3f e\n", sum.charge);
printf("Total charge on molecule 1 interface region: %11.3f e\n", sum.charge_mol_1);
printf("Total charge on molecule 2 interface region: %11.3f e\n", sum.charge_mol_2);
printf("Sum of products of contacting atoms charges at interface: %11.3f  e^2\n", sum.ccsum);
printf("Pearson correlation of contacting atoms charges at interface: %6.3f\n", sum.ccpc);
        printf("--------------------------------\n\n");

printf("Interaction energy\n");
printf("--------------------------------\n");
printf("Total Coulomb   energy: %14.4lf (kJ/mol)\n", sum.coul_nrg);
printf("Total solvation energy: %14.4lf (kJ/mol)\n", sum.solv_nrg);
printf(" --Born self    energy: %14.4lf (kJ/mol)\n", sum.born_nrg);
printf(" --Coulomb solv energy: %14.4lf (kJ/mol)\n", sum.coul_solv_nrg);
printf(" --sasa         energy: %14.4lf (kJ/mol)\n", sum.sasa_nrg);
printf("Total           energy: %14.4lf (kJ/mol)\n", (sum.coul_nrg + sum.solv_nrg));
printf("--------------------------------\n\n");

        printf("Surface buried at interface\n");
        printf("--------------------------------\n");
printf("Solvent accessible surface area buried at interface: %11.3f A^2\n", sum.sas);
if(flag_par.msms || flag_par.ns)
printf("Solvent excluded surface area buried at interface: %11.3f A^2\n", sum.ses);
printf("--------------------------------\n\n");

if(flag_par.msms || flag_par.ns)
{
        printf("Electrostatic potential at interface\n");
        printf("--------------------------------\n");
printf("Sum of products of contacting surface points potentials at interface: %11.3f (kJ/ (e mol))^2\n", sum.ppsum);
printf("Pearson correlation of contacting surface points potentials at interface: %6.3f\n", sum.pppc);
}
else
{
        printf("No information on Electrostatic potential at interface\n");
        printf("because nor NanoShaper or msms was used to generate the solvent excluded surface\n");
}
fclose(fpsum);
}

