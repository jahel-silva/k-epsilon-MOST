#include "udf.h"
#include "sg.h"
#include "mem.h"
#include <math.h>

#define L -200.0
#define loc 2
#define z0 3.32e-06
#define ufric 0.194
#define cmu 0.033
#define se 1.3
#define ce1 1.1414
#define ce2 1.92
#define Prt 1.0
#define sk 1.0
#define zp 0.0
#define Beta 5.0
#define gamma1 16.0
#define gamma2 16.0
#define zmax 1670.0
#define CB 5.0
#define CT 0.75
#define Uref 7.5
#define u_disk 7.5
#define zref 70.0

double Phi_m(double z) 
{
    double zeta = (z + z0) / L;

    if(z/L>0)
    {
	return 1.0+Beta*zeta;
    }
    else
    {
	return pow(1.0-gamma1*zeta,-1.0/4.0);
    }	

}

double Psi_m(double z) 
{
    double zeta = (z + z0) / L;
    double phi_m = Phi_m(z);

    if(z/L>0)
    {
	return -Beta*zeta;
    }
    else
    {
	return log(1.0 / 8.0 * (1.0 + pow(phi_m,-2.0)) * (1.0 + pow(phi_m,-1.0)) * (1.0 + pow(phi_m,-1.0)))-2.0 * atan(1.0 / phi_m) + M_PI / 2.0;
    }	

}

double Phi_e(double z) 
{
    double zeta = (z + z0) / L;

    if(z/L>0)
    {
	return Phi_m(z)-zeta;
    }
    else
    {
	return 1.0-zeta;
    }	

}

double fk(double z) 
{
    if(z/L>0)
    {
	return -1.0/2.0;
    }
    else
    {
	return (-4.0 * Phi_m(z) - gamma1 * Phi_e(z) * pow(Phi_m(z),5.0)) / 8.0;
    }	

}

double GB(double z) 
{
    double zeta = (z + z0) / L;
    return - pow(ufric,3) / (KAPPA * L);
}

double fp(double z, double S2){
	double S1 = 1.0 / sqrt(cmu) * sqrt(Phi_m(z) / Phi_e(z));
	double e_value = pow(ufric,3) / (KAPPA * (z + z0)) * Phi_e(z);
	double CR = 4.5 + CB * GB(z)/e_value;
	double f0 = 1.0 + cmu * pow(S1,2) / (CR - 1.0);
	return 2.0 * f0 / (1.0 + sqrt(1.0 + 4.0 * f0 *(f0-1.0)*pow(S2 / S1,2)));
}

double Area(double a, double y1, double y2, double x1)
{
	return (pow(a,2) * asin(y2 / a) / 2.0 + y2 * sqrt(pow(a,2)-pow(y2,2))/2.0) - (pow(a,2) * asin(y1 / a) / 2.0 + y1 * sqrt(pow(a,2)-pow(y1,2))/2.0) - x1*(y2-y1);
}

double Area2(double a, double y1, double y2,double y3, double x1, double x2)
{
	return (pow(a,2) * asin(y3 / a) / 2.0 + y3 * sqrt(pow(a,2)-pow(y3,2))/2.0) - (pow(a,2) * asin(y2 / a) / 2.0 + y2 * sqrt(pow(a,2)-pow(y2,2))/2.0) - x1*(y3-y1) + x2*(y2-y1);
}

double weight(double zi, double yi) 
{
    double z = MAX(fabs(zi),1e-10);
    double y = MAX(fabs(yi),1e-10);

    double r = sqrt((z-zref)*(z-zref) + y*y);
    double theta = atan(fabs((z-zref)/y))*180/M_PI;
    double x1,x2,y1,y2,y3,w;
    if(r>45.5)
    {
	w = 0.0; 
    }
    else if(r>= 39.9)
    {
	if(theta<=45)
	{
		x1 = fabs(y)-4.0;
		y2 = sqrt(pow(40.0,2)-pow(x1,2));
		w = Area(40.0,fabs(z-zref)-4.0,MIN(y2,fabs(z-zref)+4),x1) / 64.0;	
	}
	else
	{
		x1 = fabs(z-zref)-4.0;
		y2 = sqrt(pow(40.0,2)-pow(x1,2));
		w = Area(40.0,fabs(y)-4.0,MIN(y2,fabs(y)+4),x1) / 64.0;
	}
    }
    else if(r< 35.8 & r > 35.5)
    {
	if(theta<=45)
	{
		x1 = fabs(y)-4.0;
		x2 = fabs(y)+4.0;
		y1 = fabs(z-zref)-4.0;
		y2 = sqrt(pow(40.0,2)-pow(x2,2));
		y3 = fabs(z-zref)+4.0;
		w = Area2(40.0,y1,y2,y3,x1,x2) / 64.0;	
	}
	else
	{
		x1 = fabs(z-zref)-4.0;
		x2 = fabs(z-zref)+4.0;
		y1 = fabs(y)-4.0;
		y2 = sqrt(pow(40.0,2)-pow(x2,2));
		y3 = fabs(y)+4.0;
		w = Area2(40.0,y1,y2,y3,x1,x2) / 64.0;
	}
    }
    else
    {
	w = 1.0;
    }
    return MAX(MIN(w,1),0.0);
}

DEFINE_TURBULENT_VISCOSITY(mu_t_fp,c,t)
{
	double mu_t;
	double x[ND_ND];	
	C_CENTROID(x,c,t);
	double z = x[2];
	double rho = C_R(c,t);
	double DUDX = pow(C_DUDX(c,t),2) + pow(C_DUDY(c,t),2)+ pow(C_DUDZ(c,t),2);  
	double DVDX = pow(C_DVDX(c,t),2) + pow(C_DVDY(c,t),2)+ pow(C_DVDZ(c,t),2);
	double DWDX = pow(C_DWDX(c,t),2) + pow(C_DWDY(c,t),2)+ pow(C_DWDZ(c,t),2);
	double k = C_K(c,t);
	double d = C_D(c,t);
	double S2 = k / d * sqrt(DUDX+DVDX+DWDX);
	double FP = fp(z,S2);
	C_UDMI(c,t,4) = S2;
	C_UDMI(c,t,5) = FP;

	mu_t = FP * cmu *rho * SQR(k)/d;
	return MIN(mu_t,10e5);
}

DEFINE_INIT(my_init_func,d)
{
	cell_t c;
	Thread *t;
	double xi[ND_ND],z,zeta;
	
	thread_loop_c(t,d)
	{
		begin_c_loop_all(c,t)
		{
			C_CENTROID(xi,c,t);
			z = xi[2];
			zeta = (z+z0)/L;
			C_K(c,t) = pow(ufric,2.0)/sqrt(cmu) * sqrt(Phi_e(z)/Phi_m(z));
			C_D(c,t) = pow(ufric,3.0)/(KAPPA*(z+z0)) * (Phi_e(z));
			C_U(c,t) = (ufric/KAPPA) * (log((z+z0)/z0)-Psi_m(z));
		}
		end_c_loop_all(c,t)
	}
}

DEFINE_ADJUST(UDM_storage,d)
{
	double x[ND_ND],z,zeta,phi_m, phi_e, phi_h,fe,fun,fst,CKD;
	Thread *t;
	cell_t c;

	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{	
			C_CENTROID(x,c,t);
			z = x[2];
			zeta = (z+z0)/L;
			fun = (2.0-zeta)+gamma1/2.0*(1.0-12.0*zeta+7.0*pow(zeta,2.0))-pow(gamma1,2.0)/16.0*zeta*(3.0-54.0*zeta+35.0*pow(zeta,2.0));
			fst = (2.0-zeta)-2.0*Beta*zeta*(1.0-2.0*zeta+2.0*Beta*zeta);
			CKD = SQR(KAPPA)/(sk*sqrt(cmu));
			if(z/L>0)
			{
				fe = pow(Phi_m(z),-5.0/2.0)*(2.0*Phi_m(z)-1.0);
				C_UDMI(c,t,3) = pow(ufric,3.0)/(KAPPA*L)*(-(CKD/4.0)*pow(Phi_m(z),-7.0/2.0)*pow(Phi_e(z),-3.0/2.0)*fst);
			}
			else
			{
				fe = pow(Phi_m(z),5.0/2.0)*(1.0-3.0/4.0*gamma1*zeta);
				C_UDMI(c,t,3) = pow(ufric,3.0)/(KAPPA*L)*((1.0 / zeta)*(Phi_m(z)-Phi_e(z))-1.0-CKD/4.0*pow(Phi_m(z),13.0/2.0)*pow(Phi_e(z),-3.0/2.0)*fun);
			}
			
			C_UDMI(c,t,0) = 1.0/zeta*(ce1*Phi_m(z)-ce2*Phi_e(z)+(ce2-ce1)*pow(Phi_e(z),-1.0/2.0)*fe);
			C_UDMI(c,t,1) = C_R(c,t)*GB(z);
			C_UDMI(c,t,2) = C_D(c,t)/C_K(c,t)*C_UDMI(c,t,0)*C_UDMI(c,t,1);
		}
		end_c_loop(c,t)
	}

}

DEFINE_SOURCE(gbe_source,c,t,dS,eqn)
{
	double x[ND_ND],z;
	C_CENTROID(x,c,t);
	z = x[2];
	if(z>zp)
	{
	dS[eqn] = 0.0;	
	return C_UDMI(c,t,2);
	}
	else
	{
		dS[eqn] = 0.0;	
		return 0.0; 
	} 
}

DEFINE_SOURCE(gbk_source,c,t,dS,eqn)
{
	double x[ND_ND],z;
	C_CENTROID(x,c,t);
	z = x[2];
	if(z>zp)
	{
		dS[eqn] = 0.0;
		return C_UDMI(c,t,1)-C_R(c,t)*C_UDMI(c,t,3); 
	}
	else
	{
		dS[eqn] = 0.0;	
		return 0.0; 
	}
}

DEFINE_SOURCE(thrust_source,c,t,dS,eqn)
{
	double xi[ND_ND];
	C_CENTROID(xi,c,t);
	double x = xi[0];
	double y = xi[1];
	double z = xi[2];
	double ax = 0.5 * (1.0 - sqrt(1.0 - CT));
	double UH = u_disk;
	double Acell = 8.0 * 8.0;
	double Vcell = 8.0 * 8.0 * 8.0;
	C_UDMI(c,t,6) = weight(z,y);

	if(x<4 && x>-4)
	{
		dS[eqn] = 0.0;
		return -2.0 * C_R(c,t) * pow(UH,2) * ax * (1.0 - ax) * Acell / Vcell * weight(z,y);	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;		
	}
}

DEFINE_PROFILE(x_vel,t,i)
{
    double x[ND_ND],z;
    face_t f;	

    begin_f_loop(f,t)
    {
	F_CENTROID(x,f,t);
        z = x[2];
        F_PROFILE(f,t,i) = (ufric/KAPPA) * (log((z+z0)/z0)-Psi_m(z));
    }
    end_f_loop(f,t)
}

DEFINE_PROFILE(k_prof,t,i)
{
    double x[ND_ND],z;
    face_t f;	

    begin_f_loop(f,t)
    {
	F_CENTROID(x,f,t);
        z = x[2];
        F_PROFILE(f,t,i) = pow(ufric,2.0)/sqrt(cmu) * sqrt(Phi_e(z)/Phi_m(z));
    }
    end_f_loop(f,t)
}

DEFINE_PROFILE(e_prof,t,i)
{
    double x[ND_ND],z;
    face_t f;	

    begin_f_loop(f,t)
    {
	F_CENTROID(x,f,t);
        z = x[2];
        F_PROFILE(f,t,i) = pow(ufric,3.0)/(KAPPA*(z+z0)) * (Phi_e(z));
    }
    end_f_loop(f,t)
}

DEFINE_WALL_FUNCTIONS(user_log_law, f, t, c0, t0, wf_ret, yPlus, Emod)
{
	double zplus,x[ND_ND],z,C_NU,Eprime,ustar_ground;
	double wf_value;
	C_CENTROID(x,c0,t0);
	z = x[2];
	C_NU = C_MU_L(c0,t0) / C_R(c0,t0);
	ustar_ground = pow(C_K(c0,t0),1.0/2.0)*pow(cmu, 1.0/4.0);
	zplus = (z + z0) * ustar_ground / C_NU;
	Eprime = C_NU/(z0*ustar_ground);
	
	switch (wf_ret)
	{
	case UPLUS_LAM:
		wf_value = zplus;
		break;
	case UPLUS_TRB:
		wf_value = 1.0 / KAPPA * log(Eprime*zplus);
		break;
	case DUPLUS_LAM:
		wf_value = 1.0;
		break;
	case DUPLUS_TRB:
		wf_value = 1.0 / KAPPA * (1.0 / zplus);
		break;
	case D2UPLUS_TRB:
		wf_value = -1.0 / KAPPA * SQR(1.0 / zplus);
		break;
	default:
		printf("Wall function return value unavailable\n");
	}
	return wf_value;
}

DEFINE_PROFILE(zs,t,i)
{
    double x[ND_ND],z,E,Cs;
    E = 9.793;
    Cs = 0.5;
    face_t f;	

    begin_f_loop(f,t)
    {
        F_PROFILE(f,t,i) = E/Cs*z0;
    }
    end_f_loop(f,t)
}

DEFINE_SOURCE(e_dirichlet_top,c,t,dS,eqn)
{
	double x[ND_ND],z,efix,Su,Sp;
	C_CENTROID(x,c,t); 
	z = x[2];
	efix = pow(ufric,3.0) / (KAPPA*(zmax+z0))*Phi_e(zmax);
	Su = 10e5*efix;
	Sp = -10e5;
	if(z>zmax)
	{
		dS[eqn] = Sp;
		return Sp*C_D(c,t)+Su;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}

DEFINE_SOURCE(u_dirichlet_top,c,t,dS,eqn)
{
	double x[ND_ND],z,ufix,Su,Sp;
	C_CENTROID(x,c,t); 
	z = x[2];
	ufix = ufric/KAPPA * (log((zmax+z0)/z0)-Psi_m(zmax));
	Su = 10e5*ufix;
	Sp = -10e5;
	if(z>zmax)
	{
		dS[eqn] = Sp;
		return Sp*C_U(c,t)+Su;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}

DEFINE_SOURCE(k_dirichlet_top,c,t,dS,eqn)
{
	double x[ND_ND],z,kfix,Su,Sp;
	C_CENTROID(x,c,t); 
	z = x[2];
	kfix = pow(ufric,2.0) / sqrt(cmu)*sqrt(Phi_e(zmax)/Phi_m(zmax));
	Su = 10e5*kfix;
	Sp = -10e5;
	if(z>zmax)
	{
		dS[eqn] = Sp;
		return Sp*C_K(c,t)+Su;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}

DEFINE_SOURCE(e_neumann_top,c,t,dS,eqn)
{
	double x[ND_ND],z,DEDZ,Acell,Vcell;
	C_CENTROID(x,c,t); 
	z = x[2];
	Acell = 10 * 10;
	Vcell = 10 * 10 * 89.313;
	DEDZ = -pow(ufric,3) / (KAPPA * pow(z+z0,2));
	if(z>zmax)
	{
		dS[eqn] = 0.0;
		return DEDZ * C_R(c,t) * ufric * KAPPA * (z + z0) / (se * Phi_m(z)) * Acell / Vcell;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}

DEFINE_SOURCE(u_neumann_top,c,t,dS,eqn)
{
	double x[ND_ND],z,Acell,Vcell,DUDZ;
	C_CENTROID(x,c,t); 
	z = x[2];
	Acell = 10 * 10;
	Vcell = 10 * 10 * 89.313;
	DUDZ = ufric / (KAPPA * (z+z0)) * Phi_m(z);
	if(z>zmax)
	{
		dS[eqn] = 0.0;
		return DUDZ * C_R(c,t) * ufric * KAPPA * (z + z0) / Phi_m(z) * Acell / Vcell;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}

DEFINE_SOURCE(k_neumann_top,c,t,dS,eqn)
{
	double x[ND_ND],z,Acell,Vcell,DKDZ;
	C_CENTROID(x,c,t); 
	z = x[2];
	Acell = 10 * 10;
	Vcell = 10 * 10 * 89.313;
	if(z>zmax)
	{
		dS[eqn] = 0.0;
		return  C_R(c,t) * pow(ufric,3.0) * KAPPA / sqrt(cmu) * (z + z0) / L * 1.0 / pow(Phi_m(z),3.0) * pow(Phi_e(z)/Phi_m(z),-1.0/2.0) * fk(z) * Acell / Vcell;	
	}
	else
	{
		dS[eqn] = 0.0;
		return 0.0;	
	}
	
}