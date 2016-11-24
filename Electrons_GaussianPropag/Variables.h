//*****************
// LENGTH IS IN UM
// TIME IS IN FS
//*****************
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>
#include <string.h>
#include <omp.h>

using namespace std;

#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED

#define c  0.299792458
#define pi 3.141592654

//*****************
// 3D VECTOR CLASS
//*****************
class Vector
{
  public:
    double x,y,z;
    Vector () { x=0; y=0; z=0; };
    Vector (double x0, double y0, double z0) { x=x0; y=y0; z=z0; };
    Vector operator + (Vector u){ return Vector(x+u.x,y+u.y,z+u.z); }
    Vector operator - (Vector u){ return Vector(x-u.x,y-u.y,z-u.z); }
    Vector operator * (Vector u){ return Vector(x*u.x,y*u.y,z*u.z); }
    Vector operator * (double u){ return Vector(x*u,y*u,z*u); }
   	Vector operator / (double u){ return Vector(x/u,y/u,z/u); }
    Vector Cross(Vector u) { return Vector(y*u.z-z*u.y,z*u.x-x*u.z,x*u.y-y*u.x); }
    double Scalar(Vector u)	{ return (x*u.x+y*u.y+z*u.z); }
    double Length() 		{ return sqrt(x*x+y*y+z*z); }
};

//****************
// EM FIELD CLASS
//****************
class EMfield {
	public:
		EMfield() { };
		EMfield(Vector U, Vector V) { E=U; B=V; };
        EMfield operator + (EMfield EM) { return EMfield( E+EM.E, B+EM.B ) ; };
        EMfield operator - (EMfield EM) { return EMfield( E-EM.E, B-EM.B ) ; };
    public:
		Vector E,B;
};


//*******************
// LASER PULSE CLASS
//*******************
class LaserPulse {
  public:
  	LaserPulse(double a0, double w0, double ell, double t0, double phi, double lambda, double vp, Vector pol, Vector P, int fo)
  	:m_a0(a0), m_w0(w0), m_ell(ell), m_t0(t0), m_phi(phi), m_lambda(lambda), m_vp(vp), e1(pol/pol.Length()), e3(P/P.Length()), m_fo(fo)
  	{
			m_f = 2*pi*c/m_lambda; 				// frequency
			m_k = 2*pi/m_lambda/m_vp;			// wave number
  		    m_E0 = m_a0*m_f;					// Electric field;
            m_w0x = m_w0;
            m_w0y = m_w0 * ell;
			m_zrx = pi*m_w0x*m_w0x/m_lambda; // Rayleigh length
			m_zry = pi*m_w0y*m_w0y/m_lambda; // Rayleigh length

			e2 = e3.Cross(e1);
			if (abs(e1.Scalar(e3)) > 1e-6)
				cout << "ERROR: Polarization not perpendicular to propagation!" << abs(e1.Scalar(e3))<<endl;
//			eps  = (0<m_fo) * m_w0/m_zr;
            epsx = (0<m_fo) * m_w0x/m_zrx; // M.T.
            epsy = (0<m_fo) * m_w0y/m_zry; // M.T.
            eps  = (0<m_fo) * sqrt(epsx * epsy); // M.T.
			eps2 = (1<m_fo) * eps*eps;
			eps3 = (2<m_fo) * eps*eps2;
			eps4 = (3<m_fo) * eps*eps3;
			eps5 = (4<m_fo) * eps*eps4;
  	};

 // M.T.
    Vector GetW(Vector R, double t){
      double wx;
      double wy;
      Vector V;
			wx = m_w0x * sqrt( 1 + pow( R.z/m_zrx, 2) );
			wy = m_w0y * sqrt( 1 + pow( R.z/m_zry, 2) );
      V.x = wx;
      V.y = wy;
      return V;
    }



// M.T.
    double GetEnvelope(Vector R, double t){
			double x,y,z,r,w;
			double xi,xi2,nu;
			double rho,rho2,rho4,rho6,rho8;
			double phi0,phip,phir,phig,phi;
			double s,s2,s3,s4,s5,s6,s7;
			double S,S2,S3,S4,S5,S6;
			double C,C2,C3,C4,C5,C6,C7;
			double E,E1,E2,E3,B1,B2,B3;
      double wx, wy, phirx, phiry, phigx, phigy; // M.T.
      double sx, sy, Cx, Cy; // M.T.

			x = e1.Scalar(R);
			y = e2.Scalar(R);
			z = e3.Scalar(R);
			r = sqrt(x*x+y*y);

			wx = m_w0x * sqrt( 1 + pow( z/m_zrx, 2) ); // M.T.
			wy = m_w0y * sqrt( 1 + pow( z/m_zry, 2) ); // M.T.

			phi0 = pi/2+m_phi;
			phip = m_f*t - m_k*z;
      phirx = (z==0) ? 0 : x*x*z/(wx*wx*m_zrx); // M.T.
      phiry = (z==0) ? 0 : y*y*z/(wy*wy*m_zry); // M.T.
//      phirx = (z==0) ? 0 : m_k/2*x*x/m_zrx/m_zrx/(1+z*z/m_zrx/m_zrx); // M.T.
//      phiry = (z==0) ? 0 : m_k/2*y*y/m_zry/m_zry/(1+z*z/m_zry/m_zry); // M.T.

			phigx = atan(z/m_zrx); // M.T.
			phigy = atan(z/m_zry); // M.T.
			phi = phi0 + phip - phirx - phiry + phigx/2 + phigy/2;

			xi = x/m_w0x; // M.T.
			nu = y/m_w0y; // M.T.
			sx = m_w0x/wx; // M.T.
			sy = m_w0y/wy; // M.T.
			S =		sin(phi);
			Cx = sx*cos(phi+phigx); // M.T.
			Cy = sy*cos(phi+phigy); // M.T.

			E = sqrt(sx*sy)*exp(-pow(x/wx,2))*exp(-pow(y/wy,2)); // M.T.
//			E = sqrt(sx*sy)*exp(-pow(x/wx,2))*exp(-pow(y/wy,2))*exp(-pow((t-z/c/m_vp)/m_t0,2)); // M.T.

//      double Env;

//			Env = sqrt(m_w0x/m_w0x * sqrt( 1 + pow( e3.Scalar(R)/m_zrx, 2) )*m_w0y/m_w0y * sqrt( 1 + pow( e3.Scalar(R)/m_zry, 2) ))*exp(-pow(e1.Scalar(R)/m_w0x * sqrt( 1 + pow( e3.Scalar(R)/m_zrx, 2) ),2))*exp(-pow(e2.Scalar(R)/m_w0y * sqrt( 1 + pow( e3.Scalar(R)/m_zry, 2) ),2))*exp(-pow((t-e3.Scalar(R)/c/m_vp)/m_t0,2)); // M.T.

      return E;
    }


		EMfield GetEMfield(Vector R, double t)
		{
			double x,y,z,r,w;
			double xi,xi2,nu;
			double rho,rho2,rho4,rho6,rho8;
			double phi0,phip,phir,phig,phi;
			double s,s2,s3,s4,s5,s6,s7;
			double S,S2,S3,S4,S5,S6;
			double C,C2,C3,C4,C5,C6,C7;
			double E,E1,E2,E3,B1,B2,B3;
      double wx, wy, phirx, phiry, phigx, phigy; // M.T.
      double sx, sy, Cx, Cy; // M.T.

			x = e1.Scalar(R);
			y = e2.Scalar(R);
			z = e3.Scalar(R);
			r = sqrt(x*x+y*y);
            // <0 focus in electron beam >0: focus before
			//z = z + c*(50);

			if(m_fo<=1) {

			wx = m_w0x * sqrt( 1 + pow( z/m_zrx, 2) ); // M.T.
			wy = m_w0y * sqrt( 1 + pow( z/m_zry, 2) ); // M.T.

			phi0 = pi/2+m_phi;
			phip = m_f*t - m_k*z;
      phirx = (z==0) ? 0 : x*x*z/(wx*wx*m_zrx); // M.T.
      phiry = (z==0) ? 0 : y*y*z/(wy*wy*m_zry); // M.T.
//      phirx = (z==0) ? 0 : m_k/2*x*x/m_zrx/m_zrx/(1+z*z/m_zrx/m_zrx); // M.T.
//      phiry = (z==0) ? 0 : m_k/2*y*y/m_zry/m_zry/(1+z*z/m_zry/m_zry); // M.T.

			phigx = atan(z/m_zrx); // M.T.
			phigy = atan(z/m_zry); // M.T.
			phi = phi0 + phip - phirx - phiry + phigx/2 + phigy/2;

			xi = x/m_w0x; // M.T.
			nu = y/m_w0y; // M.T.
			sx = m_w0x/wx; // M.T.
			sy = m_w0y/wy; // M.T.
			S =		sin(phi);
			Cx = sx*cos(phi+phigx); // M.T.
			Cy = sy*cos(phi+phigy); // M.T.[-Npz/2:Npz/2]

// cout << phi << " " << phi0 << " " << phip << " " << phigx << " " << phigy << endl;

			E = m_E0*sqrt(sx*sy)*exp(-pow(x/wx,2))*exp(-pow(y/wy,2))*exp(-pow((t-z/c/m_vp)/m_t0,2)); // M.T.

//    printf("%f\n",m_E0);

				E1 = E*S;
				E2 = 0;
				E3 = E*xi*eps*Cx;
				B1 = 0;
				B2 = E*S;
				B3 = E*nu*eps*Cy;
			} else {
			w = m_w0 * sqrt( 1 + pow( z/m_zr, 2) );

			phi0 = pi/2+m_phi;
			phip = m_f*t - m_k*z;
			phir = (z==0) ? 0 : m_k*r*r/2/(z+m_zr*m_zr/z);

			phig = atan(z/m_zr);
			phi = phi0 + phip - phir + phig;

			xi = x/m_w0;
			nu = y/m_w0;
			s = m_w0/w;
			S =		sin(phi);
			C = s*cos(phi+phig);

			E = m_E0*s*exp(-pow(r/w,2))*exp(-pow((t-z/c/m_vp)/m_t0,2));

				xi2 = xi*xi;
				rho  = r/m_w0;
				rho2 = rho*rho;
				rho4 = rho2*rho2;
				s2 = s*s;
				s3 = s*s2;
				s4 = s*s3;
				S2 = s2*sin(phi + 2*phig);
				S3 = s3*sin(phi + 3*phig);
				C2 = s2*cos(phi + 2*phig);
				C3 = s3*cos(phi + 3*phig);
				C4 = s4*cos(phi + 4*phig);
			if(m_fo<=3) {

				E1 = E*(S + eps2*(xi2*S2-rho4*S3/4) );
				E2 = E*xi*nu*eps2*S2;
				E3 = E*xi*(eps*C + eps3*( -C2/2 + rho2*C3 - rho4*C4/4));
				B1 = 0;
				B2 = E*(S + eps2*(rho2*S2/2 - rho4*S3/4));
				B3 = E*nu*( eps*C + eps3*(C2/2 + rho2*C3/2 - rho4*C4/4) );
			} else {
				rho6 = rho2*rho4;
				rho8 = rho2*rho6;
				s5 = s*s4;
				s6 = s*s5;
				s7 = s*s6;
				S4 = s4*sin(phi + 4*phig);
				S5 = s5*sin(phi + 5*phig);
				S6 = s6*sin(phi + 6*phig);
				C5 = s5*cos(phi + 5*phig);
				C6 = s6*cos(phi + 6*phig);
				C7 = s7*cos(phi + 7*phig);
			if(m_fo<=5) {
				E1 = E*(S + eps2*(xi2*S2-rho4*S3/4) + eps4*(S2/8 - rho2*S3/4 - (rho4-16*rho2*xi2)*S4/16 - (rho6+2*xi2*rho4)*S5/8 + rho8*S6/32) );
				E2 = E*xi*nu*(eps2*S2 + eps4*(rho2*S4-rho4*S5/4) );
				E3 = E*xi*(eps*C + eps3*( -C2/2 + rho2*C3 - rho4*C4/4) + eps5*( -3*C3/8 - 3*rho2*C4/8 + 17*rho4*C5/16 - 3*rho6*C6/8 + rho8*C7/32 ) );
				B1 = 0;
				B2 = E*(S + eps2*(rho2*S2/2 - rho4*S3/4) + eps4*(-S2/8 + rho2*S3/4 + 5*rho4*S4/16 - rho6*S5/4 + rho8*S6/32) );
				B3 = E*nu*( eps*C + eps3*(C2/2 + rho2*C3/2 - rho4*C4/4) + eps5*(3*C3/8 + 3*rho2*C4/8 + 3*rho4*C5/16 - rho6*C6/4 + rho8*C7/32) );
			} } }
			if(R.z>=0){
			return EMfield(e1*E1 + e2*E2 + e3*E3, e1*B1 + e2*B2 + e3*B3);
			}
			else return EMfield(e1*0 + e2*0 + e3*0, e1*0 + e2*0 + e3*0);
		}
  private:
		double m_a0;			// Normalized vector potential
    double m_w0; 			// 1/e^2 radius
    double m_ell;
    double m_t0; 			// 1/e^2 half duration
		double m_E0; 			// Peak E-field (normalized to m_E0=qE/mc)
		double m_phi;			// phase
  	double m_lambda;	// wavelength
		double m_k;				// wavenumber
		double m_f;				// angular frequency
		double m_zr; 			// Rayleigh length
		double m_vp;			// Phase velocity
		Vector e1;				// Polarization vector
		Vector e3;				// Propagation vector
		int 	 m_fo;			// Field order

    double epsx, epsy, m_w0x, m_w0y, m_zrx, m_zry;

	private:
		Vector e2;				// Magnetic vector
		double eps,eps2,eps3,eps4,eps5;
};


//************
// WAKE CLASS
//************
class Wake {
	public:
		Wake(const char * wkdir, double v, Vector P)
		:vp(v), ez(P/P.Length())
		{
			ifstream fr,fz,fEr,fEz,fBy;
			bool err;
			char fname[100];
			unsigned int idz,idr;
			double a;
			fr.open(strcat(strcpy(fname,wkdir),"r"),ios::in);
			fz.open(strcat(strcpy(fname,wkdir),"z"),ios::in);
			fEr.open(strcat(strcpy(fname,wkdir),"Er"),ios::in);
			fEz.open(strcat(strcpy(fname,wkdir),"Ez"),ios::in);
			fBy.open(strcat(strcpy(fname,wkdir),"By"),ios::in);
			err=fr.fail()||fz.fail()||fEr.fail()||fEz.fail()||fBy.fail();
			if(!err) {
				for(fr >> a; !fr.eof(); fr >> a) r.push_back(a);
				for(fz >> a; !fz.eof(); fz >> a) z.push_back(a);
				Dr=(r.back()-r.front())/(r.size()-1);
				Dz=(z.back()-z.front())/(z.size()-1);
				Er.resize(r.size());
				Ez.resize(r.size());
				By.resize(r.size());
				for(idr=0;idr<r.size();idr++) {
					Er[idr].resize(z.size());
					Ez[idr].resize(z.size());
					By[idr].resize(z.size());
					for(idz=0;idz<z.size();idz++) {
						fEr >> Er[idr][idz];
						fEz >> Ez[idr][idz];
						fBy >> By[idr][idz];
					}
			 	}
				fr.close(); fz.close(); fEr.close(); fEz.close(); fBy.close();
		 	} else cout<<"FAILED TO OPEN WAKE DATA IN "<<wkdir<<endl;
		};
		EMfield GetEMfield(Vector R,double t) {bool keep; return GetEMfield(R,t,keep);};

		EMfield GetEMfield(Vector R,double t, bool& keep)
		{
			EMfield EM;
			int idr,idz;
			double ri,zi,dr,dz;
			double Eri,Ezi,Byi;
			Vector er,ey;

			zi = ez.Scalar(R);
			er = R - ez*zi;
			zi = zi - c*vp*t;
			ri = er.Length();
			if(ri!=0) er=er/ri;
			ey = ez.Cross(er);
			if(zi>=z.front() && ri<=r.back() )	keep=true;
			if(ri<r.front() || zi<z.front() || ri>=r.back() || zi>=z.back()) return EM;
			idr = (ri-r.front())/Dr;
			idz = (zi-z.front())/Dz;
		  dr=(ri-r[idr])/Dr;
		  dz=(zi-z[idz])/Dz;
			Eri=Er[idr][idz]*(1.-dr)*(1.-dz)+ Er[idr+1][idz]*dr*(1.-dz)+ Er[idr][idz+1]*(1.-dr)*dz+ Er[idr+1][idz+1]*dr*dz;
			Ezi=Ez[idr][idz]*(1.-dr)*(1.-dz)+ Ez[idr+1][idz]*dr*(1.-dz)+ Ez[idr][idz+1]*(1.-dr)*dz+ Ez[idr+1][idz+1]*dr*dz;
			Byi=By[idr][idz]*(1.-dr)*(1.-dz)+ By[idr+1][idz]*dr*(1.-dz)+ By[idr][idz+1]*(1.-dr)*dz+ By[idr+1][idz+1]*dr*dz;
			EM.E = er*Eri + ez*Ezi;
			EM.B = ey*Byi;
			return EM;
		};
	private:
		double vp; // propagation velocity
		Vector ez; // propagation vector
		double Dr,Dz; // grid step
		std::vector<double> r,z; // tables for coordinates
		std::vector<std::vector<double> > Er,Ez,By; // tables for fields
};

//****************
// PARTICLE CLASS
//****************
class Particle {
	public:
		Particle() {R=Vector(0,0,0); V=Vector(0,0,0); Ep=Vector(0,0,0); Bp=Vector(0,0,0); G=1.; Q=0; M=0; isout=false; nbSignChange=0; t0=-1000;};
		Particle(Vector r, Vector v,double g, double q, double m, int id,double tin)
		:R(r), V(v), Ep(Vector()), Bp(Vector()), G(g), Q(q), M(m), ID(id), isout(false), nbSignChange(0),t0(tin)
		{}
	public:
  	Vector R;
  	Vector V;
    Vector Ep;
    Vector Bp;
  	double G;
  	double Q;
  	double M;
  	int ID;
    bool isout;
    int nbSignChange;
    // mai variables:
    double t0;
    //bool isNot; // turns to 1 when particle x<0 has been crossed
};

//***************************
// PARTICLE POPULATION CLASS
//***************************
class ParticlePopulation {
	public:
		std::list<Particle> pp;
		double rnd(double min, double max) { return( min+(max-min)*rand()/RAND_MAX ); };
	public:
	        ParticlePopulation()
	    {
	         ifstream fin;
	        bool err=false;
	        char str[256];
	        double x,y,z,vx,vy,vz,t0;
	        int k = 0;

            fin.open("PopulationFile.par");

            while(!fin.getline(str,256).eof() && !err){
           // printf("%s\n",str); // printing line
            err|=7!= sscanf(str,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&x,&y,&z,&vx,&vy,&vz,&t0);

            // update particule list:
            Particle P(Vector(x,y,z),Vector(vx,vy,vz),1,-1,1,0,t0);
            P.ID = k;
            pp.push_back(P);
            k++;
                                                      }


                            };

		ParticlePopulation(double charge, double mass, int num, Vector L, Vector R0, Vector LP, Vector P0, int SysCoord)
		{
			Particle P(Vector(0,0,0),Vector(0,0,0),1,charge,mass,0,-1000);

      double theta,r,phi,thetav,rv,u1,u2,u3,u4;
			for(int k=0; k<num; k++) {
                if (L.x <= 0. && L.y <= 0.){
                  theta = rnd(0,2*pi);
                  r = rnd(0,1);
                  r = sqrt(r)*(-L.x);
                  P.R = R0+Vector(r*cos(theta),r*sin(theta)*(-L.y),rnd(-L.z,L.z));
                }
                else if (L.x >= 0. && L.y >= 0.){
          				P.R = R0+Vector(rnd(-L.x,L.x),rnd(-L.y,L.y),rnd(-L.z,L.z));
                }
                else if (L.x >= 0. && L.y <= 0.){
                  r = rnd(0,1);
                  theta = rnd(0,1);
                  P.R = R0 + Vector(L.x*sqrt(-2*log(r)) * cos(2*pi*theta),-L.y*sqrt(-2*log(r)) * sin(2*pi*theta),rnd(-L.z,L.z));
                }
                else {
                  printf("ERREUR : L.x et L.y doivent être de même signe :\n");
                  printf(" * >0 pour coordonnées cartésiennes (x,y,z)\n");
                  printf(" * <0 pour distribution elliptique (r,ell,z)\n");
                  exit(0);
                }
                if (SysCoord==0){
                  P.V = P0 + Vector(rnd(-LP.x,LP.x),rnd(-LP.y,LP.y),rnd(-LP.z,LP.z));
                }
                if (SysCoord==1){
                  rv = P0.x + LP.x*rnd(-1,1);
                  thetav = P0.y + LP.y*rnd(-1,1);
                  P.V = Vector(rv*sin(thetav),0,rv*cos(thetav));
                }
                if (SysCoord==2){
                  u1 = rnd(0,1);
                  u2 = rnd(0,1);
                  u3 = rnd(0,1);
                  u4 = rnd(0,1);
                  P.V = P0 + Vector(sqrt(-2*log(u1))*cos(2*pi*u2)*LP.x,sqrt(-2*log(u1))*sin(2*pi*u2)*LP.y,sqrt(-2*log(u3))*cos(2*pi*u4)*LP.z);
                }
                P.ID = k;
				pp.push_back(P);
			}
		};

		void MoveParticles(std::vector<LaserPulse> lpv, std::vector<Wake> wkv, double t, double dt)
		{
		  double rqm,t_rsg,a1,a2,a3,a4,a5,a6,det;
		  double rb1g,rb2g,rb3g, T;
		  Vector T_P,T_P2, EF, BF, W;
			Particle P,P1;
	  	EMfield EM,EM0;
			int lpvsz=lpv.size();
            int wkvsz=wkv.size();
			int k;
			std::list<Particle>::iterator ppit;
	  	bool keep;

		  for (ppit=pp.begin(); ppit!=pp.end(); ppit++)	{
		  	// Get particle from population
		  	P= (*ppit);

		  	EM=EM0; keep=(wkvsz==0);
				for(k=0;k<lpvsz;k++)	EM=EM + lpv.at(k).GetEMfield(P.R,t);
				for(k=0;k<wkvsz;k++)	EM= EM + wkv.at(k).GetEMfield(P.R,t,keep);

				// if outside wkfield, remove particle from population
				if(!keep) ppit = (--pp.erase(ppit));
				else {
				// First half-action by E-field
			  rqm = 0.5*dt*P.Q/P.M;
			  EF  = EM.E * rqm;

        P1.Ep = EM.E;
        P1.Bp = EM.B;
        P1.t0 = P.t0;

        T_P = P.V + EF;
			  // Magnetic rotation
			  t_rsg = rqm / sqrt(1.+T_P.Scalar(T_P)); // rqm/gamma at instance n
				BF = EM.B * t_rsg;
			  rb1g = BF.x;
			  rb2g = BF.y;
			  rb3g = BF.z;
			  a1 = rb1g * rb1g;
		  	a2 = rb2g * rb2g;
		  	a3 = rb3g * rb3g;
		  	a4 = rb1g * rb2g;
		  	a5 = rb2g * rb3g;
		  	a6 = rb3g * rb1g;
		  	det = 1./(1.+a1+a2+a3);
		  	T_P2.x = ((1.+a1-a2-a3)*T_P.x+2.*((a4+rb3g)*T_P.y+(a6-rb2g)*T_P.z))*det;
		  	T_P2.y = ((1.+a2-a3-a1)*T_P.y+2.*((a5+rb1g)*T_P.z+(a4-rb3g)*T_P.x))*det;
		  	T_P2.z = ((1.+a3-a1-a2)*T_P.z+2.*((a6+rb2g)*T_P.x+(a5-rb1g)*T_P.y))*det;
				// Second half-action by E-field
		  	T_P = T_P2+EF;
			  // Momentum and position at n+1/2
		  	P1.G = sqrt(1.+T_P.Scalar(T_P));
		  	P1.V = T_P;
			  P1.R = P.R + P1.V/P1.G*c*dt;
		  	P1.M = P.M;
		  	P1.Q = P.Q;
        P1.ID=P.ID;

        // Tests if particle out of pulse
        W = lpv.at(k).GetW(P.R,t);
        T = lpv.at(k).GetEnvelope(P.R,t);

        if (P1.R.x*P1.R.x/(W.x*W.x) + P1.R.y*P1.R.y/(W.y*W.y) > 1 || P.isout == true){
          P1.isout = true;
        }
        else P1.isout = false;

        // Increments nbSignChange if electric field changes sign
        if ( ( (P.Ep.x < 0 & P1.Ep.x > 0) || (P.Ep.x > 0 & P1.Ep.x < 0) ) && (P1.isout == false) ) {
          P1.nbSignChange = P.nbSignChange+1;
                                                                                                   }
        else P1.nbSignChange = P.nbSignChange;

            // Update population
            if(t>=P.t0){
                  *ppit = P1;}


}
			}
		}
};


//*********************
// DATA SAVE FUNCTIONS
//*********************
void saveppv(std::vector<ParticlePopulation> ppv, int it)	{
	ofstream fout;
	char str[100];
	Particle P;
	unsigned int i;
	std::list<Particle>::iterator ppit;

	for(i=0;i<ppv.size();i++) {
		sprintf(str,"results\\part.%02d.%04d",i,it);
		fout.open(str,ios::out);
		for (ppit=ppv.at(i).pp.begin(); ppit!=ppv.at(i).pp.end(); ppit++) {
  		P=(*ppit);
			sprintf(str,"%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e",\
        P.R.x,P.R.y,P.R.z,P.V.x,P.V.y,P.V.z,P.G,P.t0);
			fout << str << endl;
		}
    fout.close();

		sprintf(str,"results\\fiel.%02d.%04d",i,it);
		fout.open(str,ios::out);
		for (ppit=ppv.at(i).pp.begin(); ppit!=ppv.at(i).pp.end(); ppit++) {
  		P=(*ppit);
			sprintf(str,"%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%i",\
        P.Ep.x,P.Ep.y,P.Ep.z,P.Bp.x,P.Bp.y,P.Bp.z,P.nbSignChange);
			fout << str << endl;
		}
    fout.close();

	}
}

void savetime(double time,int it) {
	ofstream fout;
	char str[100];
	if(it==0)	fout.open("results\\t",ios::out);
	if(it>0)	fout.open("results\\t",ios::app);
	sprintf(str,"%.5e",time);
	fout << str << endl;
	fout.close();
}

//*************************
// READ .PAR FILE FUNCTION
//*************************
bool readparam(std::vector<LaserPulse> &lpv, std::vector<Wake> &wkv,	std::vector<ParticlePopulation> &ppv, double &tmin,double &tmax,double &dt,int &save)
{
	ifstream fin;
	int num,fo=1, SysCoord;
	double eps,a0,w0,ell,t0,phi,lambda,vp=1,charge,mass;
  Vector pol,lpP,wkP,L,R0,LP,P0;
	bool err=false;
	char wkdir[256], str[256], var[3];

	LaserPulse * lp;
	Wake * wk;
	ParticlePopulation * pp;

	fin.open("3dptc.par");
	while(!fin.getline(str,256).eof() && !err)	{
		strncpy(var,str,2); var[2]='\0';
		if(strncmp(var,"SS",2)==0) {
			err|=5!=sscanf(str,"SS: eps = %lf, tmin = %lf, tmax = %lf, dt = %lf, save = %d, ",
												&eps,&tmin,&tmax,&dt,&save);
			vp=sqrt(1-eps*eps);
			if(err==true) cout<<"ERROR: Can not read SS input"<<endl;
		}
		if(strncmp(var,"LP",2)==0) {
			err|=12>sscanf(str,"LP: a0 = %lf, w0 = %lf, ell = %lf, t0 = %lf, phi = %lf, lambda = %lf, pol = (%lf, %lf, %lf), lpP = (%lf, %lf, %lf), fo = %d",
										&a0,&w0,&ell,&t0,&phi,&lambda,&pol.x,&pol.y,&pol.z,&lpP.x,&lpP.y,&lpP.z,&fo);
			if(err==true) cout<<"ERROR: Can not read LP input"<<endl;
			else {			 	lp=new LaserPulse(a0,w0,ell,t0,phi,lambda,vp,pol,lpP,fo);
								 	  lpv.push_back(*lp);
								 	  printf("LP%d: a0 = %.2f, w0 = %.1f um, ell = %.1f, t0 = %.1f fs, phi = %.2f, lambda = %.2f um, pol = (%.2f, %.2f, %.2f), lpP = (%.2f, %.2f, %.2f)",
 	  											 lpv.size(),a0,w0,ell,t0,phi,lambda,pol.x,pol.y,pol.z,lpP.x,lpP.y,lpP.z);

	                  ofstream fout;
	                  char strparam[100];
                      fout.open("results\\param",ios::out);
	                  sprintf(strparam,"%f\t%f\t%f\t%f",a0,w0,ell,t0);
	                  fout << strparam << endl;
	                  fout.close();

									 	if(fo<=5) printf(", fo = %d",fo);
									 	else cout<<endl<<"ERROR: fo>5 is not supported";
									 	printf("\n");
			}
		}
		if(strncmp(var,"WK",2)==0) {
			err|=4!=sscanf(str,"WK: wkdir = %s , wkP = (%lf, %lf, %lf)",wkdir,&wkP.x,&wkP.y,&wkP.z);
			if(err==true) cout<<"ERROR: Can not read WAKE input"<<endl;
			else {				wk=new Wake(wkdir,vp,wkP);
										wkv.push_back(*wk);
										printf("WK%d: wkdir = %s , wkP = (%.2f, %.2f, %.2f)\n",wkv.size(),wkdir,wkP.x,wkP.y,wkP.z);
			}
		}
		if(strncmp(var,"PP",2)==0) {
			err|=16!=sscanf(str,"PP: charge = %lf, mass = %lf, num = %d, L = (%lf, %lf, %lf), R0 = (%lf, %lf, %lf), LP = (%lf, %lf, %lf), P0 = (%lf, %lf, %lf), SysCoord = %i",
											&charge,&mass,&num,&L.x,&L.y,&L.z,&R0.x,&R0.y,&R0.z,&LP.x,&LP.y,&LP.z,&P0.x,&P0.y,&P0.z,&SysCoord);
			if(err==true) cout<<"ERROR: Can not read POPULATION input"<<endl;
			else {				//pp=new ParticlePopulation(charge,mass,num,L,R0,LP,P0,SysCoord);
								//modified by maimouna bocoum:
								pp=new ParticlePopulation();
										ppv.push_back(*pp);
										printf("PP%d: charge = %.1f, mass = %.1f, num = %d, L(um) = (%.1f, %.1f, %.1f), R0(um) = (%.1f, %.1f, %.1f), LP(m*c) = (%.1f, %.1f, %.1f), P0(m*c) = (%.1f, %.1f, %.1f), SysCoord = %i \n",
 			  										ppv.size(),charge,mass,num,L.x,L.y,L.z,R0.x,R0.y,R0.z,LP.x,LP.y,LP.z,P0.x,P0.y,P0.z,SysCoord);
			}
		}
	}
	if(err==false) {
		printf("eps = %.3f --> bp = %f, gp = %.2f\n",eps,vp,1/sqrt(1-vp*vp));
		printf("RUN: tmin = %.1f fs, tmax = %.1f fs, dt = %.2f fs, save = %d\n---> %.0f iterations, data will be saved %.0f times, in %.2f fs intervals\n",
							tmin,tmax,dt,save,(tmax-tmin)/dt,(tmax-tmin)/dt/(double)save, dt*(double)save);
	}
	fin.close();
	return err;
}


#endif // VARIABLES_H_INCLUDED
